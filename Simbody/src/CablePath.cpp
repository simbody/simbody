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
using namespace SimTK;

//==============================================================================
//            PATH INSTANCE INFO / POS ENTRY / VEL ENTRY
//==============================================================================

PathInstanceInfo::PathInstanceInfo
    (const Array_<CableObstacle,CableObstacleIndex>& obstacles)
:   mapObstacleToSurface(obstacles.size()),
    mapSurfaceToObstacle(0),
    obstacleDisabled(obstacles.size(),false), 
    obstaclePose(obstacles.size()) 
{
    SurfaceObstacleIndex next(0);
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        const CableObstacle::Impl& obs = obstacles[ox].getImpl();
        if (obs.isDisabledByDefault())
            obstacleDisabled[ox] = true;

        if (obs.getNumCoordsPerContactPoint() > 0) { // a surface
            mapObstacleToSurface[ox] = next++;
            mapSurfaceToObstacle.push_back(ox);
            assert(mapSurfaceToObstacle.size() == next);
        }
        obstaclePose[ox]=obs.getDefaultPoseOnBody();
    }
}

std::ostream& operator<<(std::ostream& o, const PathInstanceInfo& info) {
    o << "PathInstanceInfo nObs=" << info.getNumObstacles()
        << " nSurf=" << info.getNumSurfaceObstacles() << endl;
    o << "  obs2surf: " << info.mapObstacleToSurface << endl;
    o << "  surf2obs: " << info.mapSurfaceToObstacle << endl;
    o << "  disabled="; 
    for (CableObstacleIndex ox(0); ox < info.obstacleDisabled.size(); ++ox)
        o << " " << String(info.obstacleDisabled[ox]);
    o << "\n  pose: ";
    for (CableObstacleIndex ox(0); ox < info.obstaclePose.size(); ++ox)
        o << " " << info.obstaclePose[ox].p();
    return o << endl;
}

std::ostream& operator<<(std::ostream& o, const PathPosEntry& ppe) {
    cout << "PathPosEntry: length=" << ppe.length << endl;
    cout << "mapToActive: " << ppe.mapToActive << endl;
    cout << "mapToActiveSurface: " << ppe.mapToActiveSurface << endl;
    cout << "mapToCoords: " << ppe.mapToCoords << endl;
    cout << "eIn_G: " << ppe.eIn_G << endl;
    cout << "Fu_GB: " << ppe.Fu_GB << endl;
    cout << "witnesses=" << ppe.witnesses << endl;

    cout << "geodesics lengths=";
    for (ActiveSurfaceIndex asx(0); asx < ppe.geodesics.size(); ++asx)
        cout << ppe.geodesics[asx].getLength() << " ";
    cout << endl;
    cout << "x=" << ppe.x << endl;
    cout << "err=" << ppe.err << endl;
    cout << "J dims=" << ppe.J.nrow() << " x " << ppe.J.ncol() << endl;
    return o;
}

std::ostream& operator<<(std::ostream& o, const PathVelEntry& entry) {
    return o;
}


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

void CablePath::solveForInitialCablePath(State& state) const
{   getImpl().solveForInitialCablePath(state); }

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

Real CablePath::getIntegratedCableLengthDot(const State& state) const 
{   return getImpl().getIntegratedLengthDot(state); }

void CablePath::setIntegratedCableLengthDot(State& state, Real value) const
{   getImpl().setIntegratedLengthDot(state, value); }


//==============================================================================
//                           CABLE PATH :: IMPL
//==============================================================================

//------------------------------------------------------------------------------
//                           APPLY BODY FORCES
//------------------------------------------------------------------------------
void CablePath::Impl::
applyBodyForces(const State& state, Real tension, 
                Vector_<SpatialVec>& bodyForcesInG) const
{
    if (tension <= 0) return;

    const PathPosEntry& ppe = getPosEntry(state);
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        const ActiveObstacleIndex ax = ppe.mapToActive[ox];
        if (!ax.isValid()) // exclude disabled and inactive obstacles
            continue;
        const CableObstacle::Impl& obs = obstacles[ox].getImpl();
        const MobilizedBody& B = obs.getMobilizedBody();
        const SpatialVec& Fu_GB = ppe.Fu_GB[ax];
        B.applyBodyForce(state, tension*Fu_GB, bodyForcesInG);
    }
}



//------------------------------------------------------------------------------
//                      SOLVE FOR INITIAL CABLE PATH
//------------------------------------------------------------------------------
// TODO -- this is a stub that does nothing.
void CablePath::Impl::
solveForInitialCablePath(State& state) const {
    const PathInstanceInfo& instInfo = getInstanceInfo(state);

    // These state variables are presumed to be uninitialized here.
    PathPosEntry& ppe = updPrevPosEntry(state);
    PathVelEntry& pve = updPrevVelEntry(state);
}

//------------------------------------------------------------------------------
//                           REALIZE TOPOLOGY
//------------------------------------------------------------------------------
void CablePath::Impl::
realizeTopology(State& state) {
    // Initialize instance info from defaults.
    PathInstanceInfo instInfo(obstacles);

    cout << "In realizeTopology(): " << instInfo << endl;

    // Allocate and initialize instance state variable.
    instanceInfoIx = cables->allocateDiscreteVariable(state,
        Stage::Instance, new Value<PathInstanceInfo>(instInfo));

    // Allocate continuous variable in which to integrate Ldot; this is
    // useful as a sanity check since we should have L0+integ(Ldot)=L(t) where
    // L0 is the initial length of the cable.
    Vector initz(1, Real(0));
    integratedLengthDotIx = cables->allocateZ(state, initz);

    // Allocate cache entries for position and velocity calculations.
    Value<PathPosEntry>* posEntry = new Value<PathPosEntry>();
    posEntry->upd().setNumObstacles(instInfo.getNumObstacles(),
                                    instInfo.getNumSurfaceObstacles());
    posEntryIx = cables->allocateAutoUpdateDiscreteVariable(state,
        Stage::Velocity,    // invalidates Velocity
        posEntry,           // takes over ownership of this Value
        Stage::Position);   // update depends on positions
        
    Value<PathVelEntry>* velEntry = new Value<PathVelEntry>();
    velEntry->upd().setNumObstacles(instInfo.getNumObstacles(),
                                    instInfo.getNumSurfaceObstacles());        
    velEntryIx = cables->allocateAutoUpdateDiscreteVariable(state,
        Stage::Dynamics,    // invalidates Dynamics (forces)
        velEntry,           // takes over ownership of this Value
        Stage::Velocity);   // update depends on velocities

    // Ask the System to give us some EventIds, and ask the State for some
    // slots for the witness functions.
    const int nSurfaces = instInfo.getNumSurfaceObstacles();
    mapEventIdToObstacle.clear(); // this is an std::map, not an array
    for (SurfaceObstacleIndex sox(0); sox < nSurfaces; ++sox) {
        const EventId id = cables->getSystem().getDefaultSubsystem()
                           .createEventId(cables->getMySubsystemIndex(), state);
        mapEventIdToObstacle[id] = instInfo.mapSurfaceToObstacle[sox];
    }

    // Every surface obstacle gets an event witness function that we can use
    // to watch the cable lift off of or drop down onto the obstacle.
    eventIx.invalidate();
    if (nSurfaces) eventIx = cables->allocateEventTriggersByStage
                                        (state, Stage::Position, nSurfaces);
}

//------------------------------------------------------------------------------
//                            REALIZE INSTANCE
//------------------------------------------------------------------------------
// We now know which of the obstacles are enabled. Initialize the state
// variables assuming all enabled obstacles are active.
void CablePath::Impl::
realizeInstance(const State& state) const {
    const PathInstanceInfo& instInfo = getInstanceInfo(state);
    PathPosEntry& ppe = updPosEntry(state);
    PathVelEntry& pve = updVelEntry(state);

    assert(!instInfo.obstacleDisabled.front()); // origin
    assert(!instInfo.obstacleDisabled.back());  // termination

    ActiveObstacleIndex nActive(0);
    ActiveSurfaceIndex  nActiveSurface(0);
    int                 nx = 0; // num unknowns

    Array_<Real> xInit;
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        ppe.mapToActive[ox].invalidate();
        ppe.mapToActiveSurface[ox].invalidate();
        ppe.mapToCoords[ox] = -1;

        if (instInfo.obstacleDisabled[ox]) {
            continue; // skip disabled obstacles
        }

        // This obstacle is enabled.

        ppe.mapToActive[ox] = nActive++;
        const CableObstacle::Impl& obs = obstacles[ox].getImpl();
        const int d = obs.getNumCoordsPerContactPoint();
        if (d == 0) {
            // This obstacle is a via point (or an end point).
            continue; // skip via points
        }

        // This obstacle is a surface obstacle.

        // Assuming here that any enabled surface is active.
        ppe.mapToActiveSurface[ox] = nActiveSurface++; 

        const CableObstacle::Surface::Impl& surf =
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>(obs);

        ppe.mapToCoords[ox] = nx; // first index in x slot
        nx += 2*d; // points P and Q need coords xP and xQ

        Vec3 xPhint(0), xQhint(0);
        if (surf.hasContactPointHints())
            surf.getContactPointHints(xPhint,xQhint);
        xInit.insert(xInit.end(), &xPhint[0], &xPhint[0]+d);
        xInit.insert(xInit.end(), &xQhint[0], &xQhint[0]+d);
    }

    ppe.initialize(nActive, nActiveSurface, nx);
    pve.initialize(nActive, nActiveSurface, nx);

    for (int i=0; i<nx; ++i)
        ppe.x[i] = xInit[i];

    // For disabled (TODO: inactive) obstacles, fill in the initial guess
    // at the closest-point-to-path.
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        if (!instInfo.obstacleDisabled[ox])
            continue; // skip enabled obstacles

        //TODO: for now disabled surface is treated as inactive surface
        //and needs its closest point initialized somehow.
        const CableObstacle::Impl& obs = obstacles[ox].getImpl();
        const int d = obs.getNumCoordsPerContactPoint();
        if (d == 0) continue; // just a via point
        // Obstacle is a surface.
        const SurfaceObstacleIndex sox = instInfo.mapObstacleToSurface[ox];
        const CableObstacle::Surface::Impl& surf =
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>(obs);
        Vec3 xPhint(0), xQhint(0);
        if (surf.hasContactPointHints())
            surf.getContactPointHints(xPhint,xQhint);
        ppe.closestSurfacePoint[sox] = xPhint;
    }
    cout << "INIT PE=" << ppe << endl;
    cout << "INIT VE=" << pve << endl;
}



//------------------------------------------------------------------------------
//                            REALIZE POSITION
//------------------------------------------------------------------------------
void CablePath::Impl::
realizePosition(const State& state) const {
    ensurePositionKinematicsCalculated(state);
}



//------------------------------------------------------------------------------
//                            REALIZE VELOCITY
//------------------------------------------------------------------------------
void CablePath::Impl::
realizeVelocity(const State& state) const {
    ensureVelocityKinematicsCalculated(state);
}



//------------------------------------------------------------------------------
//                            REALIZE ACCELERATION
//------------------------------------------------------------------------------
void CablePath::Impl::
realizeAcceleration(const State& state) const {
    const PathVelEntry& pve = getVelEntry(state);

    cables->updZDot(state)[integratedLengthDotIx] = pve.lengthDot;
}



//------------------------------------------------------------------------------
//                         CALC EVENT TRIGGER INFO
//------------------------------------------------------------------------------
// The integrators call this during initialization, *after* realize(Instance).
void CablePath::Impl::
calcEventTriggerInfo(const State& state, Array_<EventTriggerInfo>& info) const
{
    std::cout << "In CablePath::Impl::calcEventTriggerInfo()\n";

    for (std::map<EventId,CableObstacleIndex>::const_iterator p =
                                        mapEventIdToObstacle.begin();
         p != mapEventIdToObstacle.end(); ++p)
    {
        EventTriggerInfo einfo;
        einfo.setEventId(p->first);
        einfo.setTriggerOnFallingSignTransition(true); // OK active->inactive
        einfo.setTriggerOnRisingSignTransition(false); // TODO: debugging
        einfo.setRequiredLocalizationTimeWindow(Real(0.1)); //10% of time scale
        info.push_back(einfo);
    }
}



//------------------------------------------------------------------------------
//                               HANDLE EVENTS
//------------------------------------------------------------------------------
void CablePath::Impl::handleEvents
    (State& state, Event::Cause cause, const Array_<EventId>& eventIds,
    const HandleEventsOptions& options, HandleEventsResults& results) const
{
    std::cout << "In CablePath::Impl::handleEvents() cause="
              << Event::getCauseName(cause) << std::endl;
    std::cout << "EventIds: " << eventIds << std::endl;


    const PathInstanceInfo& instInfo = getInstanceInfo(state);
    const int nObs = instInfo.getNumObstacles();
    PathPosEntry& currPPE = updPosEntry(state);

    // Capture the current path segments before we invalidate the state.
    // TODO: won't be a problem when we're doing this in PPE rather than
    // by disabling.
    std::map<CableObstacleIndex, Vec3> r_SQp;
    std::map<CableObstacleIndex, Vec3> r_SPn;
    for (unsigned i=0; i < eventIds.size(); ++i) {
        std::map<EventId,CableObstacleIndex>::const_iterator p =
            mapEventIdToObstacle.find(eventIds[i]);
        if (p==mapEventIdToObstacle.end())
            continue;
        const CableObstacleIndex   ox  = p->second;
        if (instInfo.obstacleDisabled[ox]) {
            Vec3 rSQp, rSPn;
            findPathSegmentForObstacle(state, instInfo, currPPE, ox, 
                                       rSQp, rSPn);
            r_SQp[ox] = rSQp;
            r_SPn[ox] = rSPn;
        }
    }

    // Also capture the xP and xQ locations for currently active objects
    // so we can restore those for the objects that don't change status.
    // And save previous geodesics.
    Array_<Vec3,CableObstacleIndex> xP(nObs,Vec3(NaN));
    Array_<Vec3,CableObstacleIndex> xQ(nObs,Vec3(NaN));
    Array_<Geodesic,CableObstacleIndex> geodesics(nObs);
    for (CableObstacleIndex ox(0); ox < nObs; ++ox) {
        const int xSlot = currPPE.mapToCoords[ox];
        if (xSlot < 0) continue;
        xP[ox] = Vec3::getAs(&currPPE.x[xSlot]);
        xQ[ox] = Vec3::getAs(&currPPE.x[xSlot+3]);

        const ActiveSurfaceIndex asx = currPPE.mapToActiveSurface[ox];
        if (!asx.isValid()) continue;
        geodesics[ox] = currPPE.geodesics[asx];
    }


    // We're going to modify the *State* variable here. The state is the
    // one used as "previous" when evaluating the current path, which is what
    // triggered the event. We'll modify the cache entry to change the 
    // obstacle activation, then write the modified entry onto the state.
    PathPosEntry& prevPPE = updPrevPosEntry(state);
    std::cout << "PrevPPE in handler=\n" << prevPPE << std::endl;
    std::cout << "CurrPPE in handler=\n" << currPPE << std::endl;

    // Look through the triggered events to see if any of the event ids 
    // belong to this cable path. If so, activate or deactive the corresponding
    // surface obstacle.
    for (unsigned i=0; i < eventIds.size(); ++i) {
        std::map<EventId,CableObstacleIndex>::const_iterator p =
            mapEventIdToObstacle.find(eventIds[i]);
        if (p==mapEventIdToObstacle.end())
            continue;
        const CableObstacleIndex   ox  = p->second;
        const SurfaceObstacleIndex sox = instInfo.mapObstacleToSurface[ox];
        const CableObstacle::Surface::Impl& obs = 
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                        (getObstacleImpl(ox));
        if (instInfo.obstacleDisabled[ox]) {
            std::cout << "ENABLE OBSTACLE " << ox << std::endl;
            // Grab points before we invalidate the state.
            Vec3 rSQp = r_SQp[ox], rSPn = r_SPn[ox];
            updInstanceInfo(state).obstacleDisabled[ox] = false;
            cables->getSystem().realize(state, Stage::Instance);
            // now it's active
            const ActiveSurfaceIndex asx = currPPE.mapToActiveSurface[ox];
            const int xSlot = currPPE.mapToCoords[ox];
            assert(xSlot >= 0); // Should have had coordinates assigned
            const UnitVec3 d(rSPn-rSQp);
            const Vec3 P = currPPE.closestSurfacePoint[sox] /*- 1e-3*d*/;
            const Vec3 Q = currPPE.closestSurfacePoint[sox] /*+ 1e-3*d*/;
            xP[ox] = P;
            xQ[ox] = Q;
            const ContactGeometry& geom = obs.getContactGeometry();
            Geodesic zeroLength;
            geom.makeStraightLineGeodesic(P,Q,d, GeodesicOptions(),zeroLength);
            geodesics[ox] = zeroLength;
        } else {
            // Save the point of last contact as the initial guess for 
            // closest point tracking.
            Vec3 P_S, Q_S; // these will be almost the same point
            obs.getContactPointsOnObstacle(state,instInfo,currPPE,P_S,Q_S);
            currPPE.closestSurfacePoint[sox] = (P_S+Q_S)/2;
            currPPE.closestPathPoint[sox] = (P_S+Q_S)/2;
            std::cout << "DISABLE OBSTACLE " << ox << std::endl;
            updInstanceInfo(state).obstacleDisabled[ox] = true;
            cables->getSystem().realize(state, Stage::Instance);
        }

        // Set coordinates xP and xQ to their former values or to their
        // initial value if a surface was enabled.
        for (CableObstacleIndex ox(0); ox < nObs; ++ox) {
            const int xSlot = currPPE.mapToCoords[ox];
            if (xSlot < 0) continue;
            const Vec3& P = xP[ox];
            const Vec3& Q = xQ[ox];
            SimTK_ASSERT1_ALWAYS(!P.isNaN() && !Q.isNaN(),
                "CablePath::Impl::handleEvents(): "
                "saved point for obstacle %d was NaN.", (int)ox); 
            Vec3::updAs(&currPPE.x[xSlot])   = P;
            Vec3::updAs(&currPPE.x[xSlot+3]) = Q;;

            const ActiveSurfaceIndex asx = currPPE.mapToActiveSurface[ox];
            if (!asx.isValid()) continue;
            currPPE.geodesics[asx] = geodesics[ox];
        }
        currPPE.witnesses[sox] = 0;
        prevPPE = currPPE; // TODO
        updPrevVelEntry(state) = updVelEntry(state); // don't force computation
    }

    std::cout << "Final CurrPPE in handler=\n" << currPPE << std::endl;

    results.setExitStatus(HandleEventsResults::Succeeded);
}

//------------------------------------------------------------------------------
//                    FIND PATH SEGMENT FOR OBSTACLE
//------------------------------------------------------------------------------
// Private helper method for finding this object's path segment (the segment
// between the previous and next obstacles).
void CablePath::Impl::
findPathSegmentForObstacle
    (const State& state, const PathInstanceInfo& instInfo, 
     const PathPosEntry& ppe, CableObstacleIndex ox, 
     Vec3& Qprev_S, Vec3& Pnext_S) const
{
    const CableObstacleIndex prevOx = ppe.findPrevActiveObstacle(ox);
    const CableObstacleIndex nextOx = ppe.findNextActiveObstacle(ox);
    assert(prevOx.isValid() && nextOx.isValid());

    const CableObstacle::Surface::Impl& thisObs = 
        SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                        (getObstacleImpl(ox));
    const CableObstacle::Impl& prevObs = getObstacleImpl(prevOx);
    const CableObstacle::Impl& nextObs = getObstacleImpl(nextOx);

    // We're going to work in the body frames Bprev, Bthis, Bnext and
    // then reexpress in S when we get to the "in" and "out" vectors.

    // Abbreviate Bp=Bprev, Bn=Bnext, B=Bthis.
    const MobilizedBody& Bp = prevObs.getMobilizedBody();
    const MobilizedBody& B  = thisObs.getMobilizedBody();
    const MobilizedBody& Bn = nextObs.getMobilizedBody();

    const Rotation R_BBp = Bp.findBodyRotationInAnotherBody(state, B);
    const Rotation R_BBn = Bn.findBodyRotationInAnotherBody(state, B);

    const Rotation& R_BpSp = prevObs.getObstaclePoseOnBody(state,instInfo).R();
    const Transform& X_BS  = thisObs.getObstaclePoseOnBody(state,instInfo);
    const Rotation& R_BS   = X_BS.R();
    const Rotation& R_BnSn = nextObs.getObstaclePoseOnBody(state,instInfo).R();
    const Rotation R_SB = ~R_BS;
    const Rotation R_SSp = R_SB*R_BBp*R_BpSp;
    const Rotation R_SSn = R_SB*R_BBn*R_BnSn;

    Vec3 Pp, Qp, Pn, Qn; // points in local body frames
    prevObs.getContactStationsOnBody(state, instInfo, ppe, Pp, Qp);
    nextObs.getContactStationsOnBody(state, instInfo, ppe, Pn, Qn);

    // Find prev and next points measured from & expressed in Bthis.
    const Vec3 r_BQp = Bp.findStationLocationInAnotherBody(state, Qp, B);
    const Vec3 r_BPn = Bn.findStationLocationInAnotherBody(state, Pn, B);

    // Now shift into S frame.
    Qprev_S = ~X_BS*r_BQp;
    Pnext_S = ~X_BS*r_BPn;
}

//------------------------------------------------------------------------------
//                  ENSURE POSITION KINEMATICS CALCULATED
//------------------------------------------------------------------------------
void CablePath::Impl::
ensurePositionKinematicsCalculated(const State& state) const {
    if (cables->isDiscreteVarUpdateValueRealized(state, posEntryIx))
        return;

    const PathInstanceInfo& instInfo = getInstanceInfo(state);
    const PathPosEntry&     prevPPE  = getPrevPosEntry(state);
    PathPosEntry&           ppe      = updPosEntry(state);

    solveForPathPoints(state, instInfo, ppe);


    //TODO: combine with earlier loop
    for (SurfaceObstacleIndex sox(0); sox < instInfo.getNumSurfaceObstacles();
         ++sox)
    {
        const CableObstacleIndex ox = instInfo.mapSurfaceToObstacle[sox];
        if (ppe.mapToActiveSurface[ox].isValid())
            continue; // skip active surface obstacles

        Vec3 r_SQp, r_SPn;
        findPathSegmentForObstacle(state, instInfo, ppe, ox, r_SQp, r_SPn);

        const Vec3 r_QpPn = (r_SPn - r_SQp);
        const UnitVec3 d_QpPn(r_QpPn); // direction of line

        const CableObstacle::Surface::Impl& thisObs = 
        SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                        (getObstacleImpl(ox));
        const ContactGeometry& geo = thisObs.getContactGeometry();

        Vec3 prevClosestPt = prevPPE.closestSurfacePoint[sox];
        if (prevClosestPt.isNaN()) {
            Vec3 xPhint, xQhint;
            if (thisObs.hasContactPointHints()) {
                thisObs.getContactPointHints(xPhint, xQhint);
                prevClosestPt = xPhint;
            } else
                prevClosestPt = (r_SPn+r_SQp)/2;
        }
        Vec3 closestPointOnSurface, closestPointOnLine; Real height;
        bool succeeded = geo.trackSeparationFromLine
            ((r_SPn+r_SQp)/2, d_QpPn,
            prevClosestPt, closestPointOnSurface, closestPointOnLine, height);

        ppe.witnesses[sox] = height;
        ppe.closestSurfacePoint[sox] = closestPointOnSurface;
        ppe.closestPathPoint[sox] = closestPointOnLine;
        //std::cout << "HEIGHT=" << ppe.witnesses[sox] << std::endl;
    }

    Vector& eventTriggers = 
        cables->updEventTriggersByStage(state,Stage::Position); 

    for (SurfaceObstacleIndex i(0); i < instInfo.getNumSurfaceObstacles(); ++i)
    {
        eventTriggers[eventIx+i] = ppe.witnesses[i];
    }

    cables->markDiscreteVarUpdateValueRealized(state, posEntryIx);
}



//------------------------------------------------------------------------------
//                  ENSURE VELOCITY KINEMATICS CALCULATED
//------------------------------------------------------------------------------
void CablePath::Impl::
ensureVelocityKinematicsCalculated(const State& state) const {
    if (cables->isDiscreteVarUpdateValueRealized(state, velEntryIx))
        return;

    const PathInstanceInfo& instInfo = getInstanceInfo(state);
    const PathPosEntry& ppe = getPosEntry(state);
    PathVelEntry&       pve = updVelEntry(state);

    pve.lengthDot = 0;
    pve.unitPower = 0;

    // Solve for xdot:
    //   Calc RHS in J xdot = Kdot.
    //   Solve xdot = J\ Kdot. 

    if (ppe.x.size()) {
        findKinematicVelocityErrors(state, instInfo, ppe, pve);
        ppe.JInv.solve(pve.nerrdotK, pve.xdot);
        cout << "*** xdot=" << pve.xdot << endl;
    }

    //TODO: calc length dot
    Vec3 vprev_GQ;
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        const ActiveObstacleIndex ax = ppe.mapToActive[ox];
        if (!ax.isValid())
            continue; // skip disabled or inactive obstacles

        const CableObstacle::Impl& obs = obstacles[ox].getImpl();
        const MobilizedBody& B = obs.getMobilizedBody();
        Vec3 P_B, Q_B, v_GP, v_GQ;
        obs.getContactStationsOnBody(state, instInfo, ppe, P_B, Q_B);

        //TODO: this is wrong -- need to use the contact point velocities
        // not the material point velocities.
        v_GP = B.findStationVelocityInGround(state, P_B);
        v_GQ = B.findStationVelocityInGround(state, Q_B);

        if (ax > 0) {
            pve.lengthDot += dot(v_GP - vprev_GQ, ppe.eIn_G[ax]);
        }
        vprev_GQ = v_GQ;

        const SpatialVec& V_GB = obs.getBodyVelocity(state);
        pve.unitPower += ~ppe.Fu_GB[ax] * V_GB;
    }

    cables->markDiscreteVarUpdateValueRealized(state, velEntryIx);
}

// Numerical Jacobian calculation.

class PathError : public Differentiator::JacobianFunction {
public:
    // Caution: we're holding references.
    PathError(int nx, 
              const CablePath::Impl&  path,
              const State&            state,
              const PathInstanceInfo& instInfo,
              const PathPosEntry&     ppe, 
              Real                    accuracy)
    :   Differentiator::JacobianFunction(nx, nx), path(path),
        state(state),instInfo(instInfo),localPpe(ppe)
    {   setEstimatedAccuracy(accuracy); }

    int f(const Vector& x, Vector& fx) const override {
        localPpe.x = x;
        path.calcPathError(state, instInfo, localPpe);
        fx = localPpe.err;
        return 0;
    }

private:
    const CablePath::Impl&  path;
    const State&            state;
    const PathInstanceInfo& instInfo;
    mutable PathPosEntry    localPpe;
};

//------------------------------------------------------------------------------
//                         PROJECT ONTO SURFACE
//------------------------------------------------------------------------------
// Take whatever is in ppe.x and move any implicit surface points that are not
// on the surface downhill to the nearest surface point.
void CablePath::Impl::
projectOntoSurface(const PathInstanceInfo& instInfo, 
                   PathPosEntry& ppe) const
{
    // Skip origin and termination "obstacles" at beginning and end.
    for (CableObstacleIndex ox(1); ox < obstacles.size()-1; ++ox) {
        const ActiveSurfaceIndex asx = ppe.mapToActiveSurface[ox];
        if (!asx.isValid())
            continue; // skip via points and inactive surfaces

        const int xSlot = ppe.mapToCoords[ox];
        assert(xSlot >= 0); // Should have had coordinates assigned

        Vec3& P = Vec3::updAs(&ppe.x[xSlot]);
        Vec3& Q = Vec3::updAs(&ppe.x[xSlot+3]);

        const CableObstacle::Surface::Impl& obs = 
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                        (getObstacleImpl(ox));

        const ContactGeometry& geom = obs.getContactGeometry();
        P = geom.projectDownhillToNearestPoint(P);
        Q = geom.projectDownhillToNearestPoint(Q);
    }
}

//------------------------------------------------------------------------------
//                         SOLVE FOR PATH POINTS
//------------------------------------------------------------------------------
// The entry unit vectors have been computed for all obstacles and the lengths
// of all straight segments have been accumulated. Now solve for the unknown
// path point locations on surface obstacles, and accumulate the resulting
// geodesic lengths to complete the path length.
void CablePath::Impl::
solveForPathPoints(const State& state, const PathInstanceInfo& instInfo, 
                   PathPosEntry& ppe) const 
{
    const PathPosEntry& prevPPE = getPrevPosEntry(state);
    if (prevPPE.x.size())
        ppe.x = prevPPE.x; // start with previous solution if there is one

    projectOntoSurface(instInfo,ppe); // clean up first

    calcPathError(state,instInfo,ppe);

    if (ppe.x.size() == 0)
        return; // only via points; no iteration to do

    const Real ftol = Real(1e-12)*1000; // TODO
    const Real xtol = Real(1e-12)*1000;

    const Real estimatedPathErrorAccuracy = ftol;
    PathError pathErrorFnc(ppe.x.size(), *this, state, instInfo, ppe, 
                           estimatedPathErrorAccuracy);
    Differentiator diff(pathErrorFnc);

    Vector dx, xold, xchg;

    Real f = ppe.err.norm();
    cout << "\n***PATH solveForPathPoints@t=" << state.getTime() 
         << " err=" << f << "\n"; 

    Real fold, lam = 1, nextlam = 1;
    Real dxnormPrev = Infinity;
    int maxIter = 20;
    for (int i = 0; i < maxIter; ++i) {
        // We always need a Jacobian even if the path is already good enough
        // because we use it to solve for xdot. So we might as well do one
        // iteration.
        if (i > 0 && f <= ftol) {
            std::cout << "\n***PATH converged in " 
                << i << " iterations err=" << f << "\n\n";
            break;
        }
        //cout << "obstacle err = " << f << ", x = " << ppe.x << endl;

        //diff.calcJacobian(ppe.x, ppe.err, ppe.J, 
        //                  Differentiator::ForwardDifference);
       // Matrix Jnum = ppe.J;
        calcPathErrorJacobian(state, instInfo, ppe);

        //cout << "DIFF J=" << Jnum-ppe.J;
        //cout << "DIFF NORM=" << (Jnum-ppe.J).norm() << "\n";

        ppe.JInv.factor(ppe.J);
        //ppe.JInv.factor(Jnum);

        fold = f;
        xold = ppe.x;

        ppe.JInv.solve(ppe.err, dx);

        const Real dxnorm = std::sqrt(dx.normSqr()/ppe.x.size()); // rms
        cout << "|dx| = " << dxnorm << endl;
        if (dxnorm > Real(.99)*dxnormPrev) {
           std::cout << "\nPATH stalled in " 
                << i << " iterations err=" << f << " |dx|=" << dxnorm << "\n\n";
            break;
        }

        // backtracking
        lam = nextlam;
        while (true) {
            xchg = lam*dx;
            ppe.x = xold - xchg;
            projectOntoSurface(instInfo,ppe); // clean up first

            calcPathError(state,instInfo,ppe);
            f = ppe.err.norm();
            //cout << "step=" << lam << " obstacle err = " << f 
            //     << ", x = " << ppe.x << endl;
            if (f <= fold)
                break;
            lam = lam / 2;
        }
        //cout << "step size=" << lam << endl;

        if (lam == nextlam)
            nextlam = std::min(2*lam, Real(1));

        dxnormPrev = dxnorm;
    }
    //cout << "obstacle error = " << ppe.err << endl;
    //SimTK_ERRCHK3_ALWAYS(f <= ftol, "CablePath::solveForPathPoints()", 
    //    "At t=%g, achieved patherr=%g but tol=%g.", state.getTime(), f, ftol);
}



//------------------------------------------------------------------------------
//                      FIND KINEMATIC VELOCITY ERRORS
//------------------------------------------------------------------------------
// Given a patherr function err(K;x), the time derivative at x=x0 is 
//    errdot(K,Kdot,x0; xdot) = Derr/Dx xdot + errdotK
//            where   errdotK = Derr/DK Kdot
// Here we calculate errdotK, which is the errdot we would get if
// xdot==0, i.e., with the contact points frozen on their surfaces. This term
// depends only on the ordinary material point kinematics. We will later use
// this to solve for the value of xdot that makes errdot==0.
//
// We actually calculate nerrdotK = -errdotK here so we don't have to negate
// it later.
void CablePath::Impl::
findKinematicVelocityErrors
   (const State& state, const PathInstanceInfo& instInfo, 
    const PathPosEntry& ppe, PathVelEntry& pve) const
{
    if (ppe.x.size() == 0)
        return; // only via points; no xdots to worry about

    // Run through all the active surface obstacles calculating time 
    // derivatives of the input unit vectors eIn and eOut and passing those 
    // to the obstacles to obtain the kinematic velocity errors. 
    // TODO: accumulate the length changes of the straight-line segments
    // here.

    // The origin point is always active. Start off with that as the "previous"
    // active obstacle.
    CableObstacleIndex prevActiveOx = CableObstacleIndex(0);

    // Now go through the active obstacles and stop to process whenever one
    // of those is a surface obstacle.
    CableObstacleIndex thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);

    // Stop when we reach the termination point -- that's not a surface.
    while (thisActiveOx < obstacles.size()-1) {
        // Is this active obstacle a surface? If not keep going.
        const ActiveSurfaceIndex asx = ppe.mapToActiveSurface[thisActiveOx];
        if (!asx.isValid()) {
            prevActiveOx = thisActiveOx;
            thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);
            continue; // skip via points and inactive surfaces
        }

        // This active obstacle is a surface. Find the next active obstacle.
        const CableObstacleIndex nextActiveOx = 
            ppe.findNextActiveObstacle(thisActiveOx);

        // Now we have prev, this, and next with prev and next active and
        // this an active surface. We want all the points in the S frame
        // and the velocities of Qprev and Pnext in the S frame (P and Q
        // velocities are zero in S).
        const CableObstacle::Impl& prevObs = getObstacleImpl(prevActiveOx);
        const CableObstacle::Surface::Impl& thisObs = 
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                (getObstacleImpl(thisActiveOx));
        const CableObstacle::Impl& nextObs = getObstacleImpl(nextActiveOx);

        // We're going to work in the body frames Bprev, Bthis, Bnext and
        // then reexpress in S when we get to the "in" and "out" vectors.
        const Rotation& R_BS = thisObs.getObstaclePoseOnBody(state,instInfo).R();
        const Rotation R_SB = ~R_BS;

        // Abbreviate Bp=Bprev, Bn=Bnext, B=Bthis.
        const MobilizedBody& Bp = prevObs.getMobilizedBody();
        const MobilizedBody& B = thisObs.getMobilizedBody();
        const MobilizedBody& Bn = nextObs.getMobilizedBody();
        Vec3 Pp, Qp, P, Q, Pn, Qn; // points in local body frames
        prevObs.getContactStationsOnBody(state, instInfo, ppe, Pp, Qp);
        thisObs.getContactStationsOnBody(state, instInfo, ppe, P, Q);
        nextObs.getContactStationsOnBody(state, instInfo, ppe, Pn, Qn);

        // Find prev and next points measured from & expressed in Bthis.
        const Vec3 r_BQp = Bp.findStationLocationInAnotherBody(state, Qp, B);
        const Vec3 r_BPn = Bn.findStationLocationInAnotherBody(state, Pn, B);

        // Find prev and next point velocities measured & expressed in Bthis.
        const Vec3 v_BQp = Bp.findStationVelocityInAnotherBody(state, Qp, B);
        const Vec3 v_BPn = Bn.findStationVelocityInAnotherBody(state, Pn, B);

        // Get relative position vectors in B and re-express in S.
        const Vec3 rIn  = R_SB*(P - r_BQp);
        const Vec3 rOut = R_SB*(r_BPn - Q);

        // These are the time derivatives in S of rIn and rOut. P, Q, and R_SB
        // are all constants so P and Q drop out.
        const Vec3 vIn  = -(R_SB*v_BQp);
        const Vec3 vOut =   R_SB*v_BPn;

        // Now create unit vectors along rIn and rOut and their time 
        // derivatives in S.
        const Real ooNormIn  = 1/rIn.norm();    // oo == "one over"
        const Real ooNormOut = 1/rOut.norm(); 

        // Create unit vectors in S; "true" here means no need to normalize.
        const UnitVec3 eIn (ooNormIn *rIn,  true);
        const UnitVec3 eOut(ooNormOut*rOut, true); 

        const Vec3 eInDot  = ooNormIn* (vIn  - (~eIn *vIn )*eIn);
        const Vec3 eOutDot = ooNormOut*(vOut - (~eOut*vOut)*eOut);

        // Now ask the obstacle to calculate the kinematic velocity errors
        // given the in/out direction time derivatives.
        
        const int xSlot = ppe.mapToCoords[thisActiveOx];
        assert(xSlot >= 0); // Should have had coordinates assigned

        Vec6::updAs(&pve.nerrdotK[xSlot]) =
            thisObs.calcSurfaceNegKinematicVelocityError
               (ppe.geodesics[asx],             // solved geodesic
                eIn,                            // in S frame
                eInDot,
                Vec3::getAs(&ppe.x[xSlot]),     // xP
                Vec3::getAs(&ppe.x[xSlot+3]),   // xQ
                eOut,                           // in S frame
                eOutDot);
        
        prevActiveOx = thisActiveOx;
        thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);       
    }
}

//------------------------------------------------------------------------------
//                            CALC PATH ERROR
//------------------------------------------------------------------------------
// Given a set of nx coordinates for the
// path points, calculate the nx error conditions that result. Return the
// new geodesics in case anyone is interested.
void CablePath::Impl::
calcPathError(const State& state, const PathInstanceInfo& instInfo, 
              PathPosEntry& ppe) const
{
    const PathPosEntry& prevPPE = getPrevPosEntry(state);

    ppe.length = 0;

    // First pass: run through all the enabled obstacles. Update the distance
    // function for inactive obstacles. Precalculate entry directions and
    // unit forces for all active obstacles, and accumulate the lengths of
    // straight-line segments.
    Vec3 prevQB_G, prevQ_G;
    for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
        if (instInfo.obstacleDisabled[ox])
            continue;

        const ActiveObstacleIndex ax = ppe.mapToActive[ox];
        if (!ax.isValid()) {
            // TODO: update distance function
            continue;
        }

        const CableObstacle::Impl& obs = getObstacleImpl(ox);
        Vec3 PB_G, QB_G;
        obs.expressContactStationsInGround(state,instInfo,ppe,
                                            PB_G, QB_G);
        const Vec3& Bo_G = obs.getBodyTransform(state).p();
        const Vec3 P_G = Bo_G + PB_G;
        const Vec3 Q_G = Bo_G + QB_G;
        if (ax == 0) { // origin
            prevQB_G = QB_G;
            prevQ_G  = Q_G;
            // No forces applied from the "entry side" at the origin.
            ppe.Fu_GB[ax] = SpatialVec(Vec3(0));
            continue;
        }

        const Vec3 r_QP_G(P_G - prevQ_G);
        const Real lineSegLen = r_QP_G.norm();
        const UnitVec3 eIn_G(r_QP_G/lineSegLen, true);
        ppe.eIn_G[ax] = eIn_G;
        ppe.length += lineSegLen;
        const SpatialVec unitForcePrevQ(prevQB_G % eIn_G,  eIn_G);
        const SpatialVec unitForceP    (   eIn_G % PB_G,  -eIn_G);
        // Add in force component from the "exit side" of previous obstacle.
        ppe.Fu_GB[ax.prev()] += unitForcePrevQ;
        // Initialize force to the contribution from the "entry side" of
        // this obstacle.
        ppe.Fu_GB[ax] = unitForceP;

        prevQB_G = QB_G;
        prevQ_G  = Q_G;
    }

    if (ppe.x.size() == 0) {
        // No active surfaces -- this path is all lines between points so 
        // has no errors.
        return;
    }

    // Pass 2 : use above info to calculate errors from active surfaces.


    // Skip origin and termination "obstacles" at beginning and end.
    for (CableObstacleIndex ox(1); ox < obstacles.size()-1; ++ox) {
        const ActiveSurfaceIndex asx = ppe.mapToActiveSurface[ox];
        if (!asx.isValid())
            continue; // skip via points and inactive surfaces

        const int xSlot = ppe.mapToCoords[ox];
        assert(xSlot >= 0); // Should have had coordinates assigned

        const CableObstacle::Surface::Impl& obs = 
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                        (getObstacleImpl(ox));

        const Rotation& R_BS = obs.getObstaclePoseOnBody(state, instInfo).R();
        const Rotation& R_GB = obs.getBodyTransform(state).R();
        const Rotation  R_GS = R_GB*R_BS;

        const ActiveObstacleIndex ax = ppe.mapToActive[ox];
        const ActiveSurfaceIndex prevASX = prevPPE.mapToActiveSurface[ox];
        const UnitVec3 eIn_S  = ~R_GS * ppe.eIn_G[ax];
        const UnitVec3 eOut_S = ~R_GS * ppe.eIn_G[ax.next()];
        Vec6::updAs(&ppe.err[xSlot]) =
            obs.calcSurfacePathError(   
                prevASX.isValid() ? prevPPE.geodesics[prevASX] : Geodesic(),
                eIn_S,
                Vec3::getAs(&ppe.x[xSlot]),  // xP
                Vec3::getAs(&ppe.x[xSlot+3]),// xQ
                eOut_S,
                ppe.geodesics[asx] );

        const Geodesic& geod = ppe.geodesics[asx];
        const Real signP = (Real)sign(dot(eIn_S,geod.getTangentP()));
        const Real signQ = (Real)sign(dot(eOut_S,geod.getTangentQ()));
        const bool isFlipped = (signP<0 && signQ<0);
        const Real geoLength = isFlipped ? -geod.getLength() : geod.getLength();

        ppe.length += geoLength;

        const SurfaceObstacleIndex sox = instInfo.mapObstacleToSurface[ox];
        assert(sox.isValid());

        ppe.witnesses[sox] = geoLength;

        //std::cout << "WITNESS=" << ppe.witnesses[sox] << std::endl;
    }

}

//------------------------------------------------------------------------------
//                          CALC PATH ERROR JACOBIAN
//------------------------------------------------------------------------------
// Assemble the nx X nx banded Jacobian J=D patherr / Dx from per-obstacle 
// blocks. The obstacles compute the Jacobian of their own path error 
// functions, which are eHat(eIn_S, xP, xQ, eOut_S), with all arguments in the
// obstacle frame S. The blocks we need are instead the Jacobian of
//    e(xQ-1, xP, xQ, xP+1) = eHat(eIn_S(xQ-1,xP), xP, xQ, eOut_S(xQ,xP+1))
// so we need to apply the chain rule terms
//          D eIn_S    D eIn_S     D eOut_S    D eOut_S
//          --------   --------    --------    --------
//           D xQ-1      D xP        D xQ       D xP+1
// to produce the block we need for the patherr Jacobian.
void CablePath::Impl::
calcPathErrorJacobian(const State&            state, 
                      const PathInstanceInfo& instInfo, 
                      PathPosEntry&           ppe) // in/out
                      const
{
    const PathPosEntry& prevPPE = getPrevPosEntry(state);

    const int nx = ppe.x.size();
    ppe.J.resize(nx, nx);
    ppe.J.setToZero();
    if (nx == 0)
        return; // only via points; nothing to do

    // Run through all the active surface obstacles calculating partial 
    // derivatives of the input unit vectors eIn and eOut and applying those
    // to the Jacobian blocks returned by the obstacles.

    // The origin point is always active. Start off with that as the "previous"
    // active obstacle.
    CableObstacleIndex prevActiveOx = CableObstacleIndex(0);

    // Now go through the active obstacles and stop to process whenever one
    // of those is a surface obstacle.
    CableObstacleIndex thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);

    // Stop when we reach the termination point -- that's not a surface.
    while (thisActiveOx < obstacles.size()-1) {
        // Is this active obstacle a surface? If not keep going.
        const ActiveSurfaceIndex asx = ppe.mapToActiveSurface[thisActiveOx];
        if (!asx.isValid()) {
            prevActiveOx = thisActiveOx;
            thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);
            continue; // skip via points and inactive surfaces
        }

        // This active obstacle is a surface. Find the next active obstacle.
        const CableObstacleIndex nextActiveOx = 
            ppe.findNextActiveObstacle(thisActiveOx);

        // Now we have prev, this, and next with prev and next active and
        // this an active surface. We want all the points in the S frame
        // and the velocities of Qprev and Pnext in the S frame (P and Q
        // velocities are zero in S).
        const CableObstacle::Impl& prevObs = getObstacleImpl(prevActiveOx);
        const CableObstacle::Surface::Impl& thisObs = 
            SimTK_DYNAMIC_CAST_DEBUG<const CableObstacle::Surface::Impl&>
                                                (getObstacleImpl(thisActiveOx));
        const CableObstacle::Impl& nextObs = getObstacleImpl(nextActiveOx);

        // We're going to work in the body frames Bprev, Bthis, Bnext and
        // then reexpress in S when we get to the "in" and "out" vectors.

        // Abbreviate Bp=Bprev, Bn=Bnext, B=Bthis.
        const MobilizedBody& Bp = prevObs.getMobilizedBody();
        const MobilizedBody& B  = thisObs.getMobilizedBody();
        const MobilizedBody& Bn = nextObs.getMobilizedBody();

        const Rotation R_BBp = Bp.findBodyRotationInAnotherBody(state, B);
        const Rotation R_BBn = Bn.findBodyRotationInAnotherBody(state, B);

        const Rotation& R_BpSp = prevObs.getObstaclePoseOnBody(state,instInfo).R();
        const Rotation& R_BS   = thisObs.getObstaclePoseOnBody(state,instInfo).R();
        const Rotation& R_BnSn = nextObs.getObstaclePoseOnBody(state,instInfo).R();
        const Rotation R_SB = ~R_BS;
        const Rotation R_SSp = R_SB*R_BBp*R_BpSp;
        const Rotation R_SSn = R_SB*R_BBn*R_BnSn;

        Vec3 Pp, Qp, P, Q, Pn, Qn; // points in local body frames
        prevObs.getContactStationsOnBody(state, instInfo, ppe, Pp, Qp);
        thisObs.getContactStationsOnBody(state, instInfo, ppe, P, Q);
        nextObs.getContactStationsOnBody(state, instInfo, ppe, Pn, Qn);

        // Find prev and next points measured from & expressed in Bthis.
        const Vec3 r_BQp = Bp.findStationLocationInAnotherBody(state, Qp, B);
        const Vec3 r_BPn = Bn.findStationLocationInAnotherBody(state, Pn, B);

        // Get relative position vectors in B and re-express in S.
        const Vec3 rIn  = R_SB*(P - r_BQp);
        const Vec3 rOut = R_SB*(r_BPn - Q);

        // Now create unit vectors along rIn and rOut, in S.
        const Real ooNormIn  = 1/rIn.norm();    // oo == "one over"
        const Real ooNormOut = 1/rOut.norm(); 

        // Create unit vectors in S; "true" here means no need to normalize.
        UnitVec3 eIn (ooNormIn *rIn,  true);
        UnitVec3 eOut(ooNormOut*rOut, true); 

        Mat33 DeinDP   =  ooNormIn *(Mat33(1) - eIn*~eIn);
        Mat33 DeoutDQ  = -ooNormOut*(Mat33(1) - eOut*~eOut);
        Mat33 DeinDQp  = -DeinDP  * R_SSp;
        Mat33 DeoutDPn = -DeoutDQ * R_SSn;

        // Now ask the obstacle to calculate the Jacobian Jhat of its error
        // function ehat(e_In, xP, xQ, e_Out).
        
        const int xSlot = ppe.mapToCoords[thisActiveOx];
        assert(xSlot >= 0); // Should have had coordinates assigned
        const ActiveSurfaceIndex prevASX = 
            prevPPE.mapToActiveSurface[thisActiveOx];
 
        //Mat63 DehatDein1, DehatDxP1, DehatDxQ1, DehatDeout1;
        //thisObs.calcSurfacePathErrorJacobianNumerically
        //   (prevASX.isValid() ? prevPPE.geodesics[prevASX] : Geodesic(),
        //    eIn,                            // in S frame
        //    Vec3::getAs(&ppe.x[xSlot]),     // xP
        //    Vec3::getAs(&ppe.x[xSlot+3]),   // xQ
        //    eOut,
        //    ppe.geodesics[asx],             // solved geodesic
        //    DehatDein1, DehatDxP1, DehatDxQ1, DehatDeout1);

        Mat63 DehatDein, DehatDxP, DehatDxQ, DehatDeout;
        thisObs.calcSurfacePathErrorJacobianAnalytically
           (eIn,                            // in S frame
            Vec3::getAs(&ppe.x[xSlot]),     // xP
            Vec3::getAs(&ppe.x[xSlot+3]),   // xQ
            eOut,
            ppe.geodesics[asx],             // solved geodesic
            DehatDein, DehatDxP, DehatDxQ, DehatDeout);

        //DehatDein=DehatDein1;DehatDxP=DehatDxP1;DehatDxQ=DehatDxQ1;DehatDeout=DehatDeout1;

        //cout << "err=" << ppe.err << "\n";

        //cout << "DehatDxP =" << DehatDxP;
        //cout << "DehatDxP1=" << DehatDxP1;
        //cout << "diff=" << (DehatDxP-DehatDxP1).norm() << ": " << (DehatDxP-DehatDxP1);

        //cout << "DehatDxQ =" << DehatDxQ;
        //cout << "DehatDxQ1=" << DehatDxQ1;
        //cout << "diff=" << (DehatDxQ-DehatDxQ1).norm() << ": " << (DehatDxQ-DehatDxQ1);

        if (xSlot >= 3)
            ppe.J(xSlot,xSlot-3,6,3) = Matrix(DehatDein*DeinDQp);
        ppe.J(xSlot,xSlot,  6,3) = Matrix(DehatDxP + DehatDein *DeinDP);
        ppe.J(xSlot,xSlot+3,6,3) = Matrix(DehatDxQ + DehatDeout*DeoutDQ);
        if (xSlot+9 < nx)
            ppe.J(xSlot,xSlot+6,6,3) = Matrix(DehatDeout*DeoutDPn);

        prevActiveOx = thisActiveOx;
        thisActiveOx = ppe.findNextActiveObstacle(prevActiveOx);       
    }

}

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
getDecorativeGeometry() const {return impl->getDecoration();}
DecorativeGeometry& CableObstacle::
updDecorativeGeometry() {return impl->updDecoration();}
bool CableObstacle::
isDisabledByDefault() const {return impl->isDisabledByDefault();}

CableObstacle& CableObstacle::
setDisabledByDefault(bool shouldBeDisabled) 
{   impl->setDisabledByDefault(shouldBeDisabled); return *this;}
CableObstacle& CableObstacle::
setDefaultTransform(const Transform& X_BS) 
{   impl->defaultX_BS=X_BS; return *this;}
CableObstacle& CableObstacle::
setDecorativeGeometry(const DecorativeGeometry& viz)
{   impl->setDecoration(viz); return *this; }

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
Surface(CablePath& path, const MobilizedBody& mobod, 
        const Transform& X_BS, const ContactGeometry& geom)
:   CableObstacle(new Surface::Impl(path, mobod, X_BS, geom))
{
    impl->index = path.updImpl().adoptObstacle(*this); 
}

CableObstacle::Surface& CableObstacle::Surface::
setNearPoint(const Vec3& stationInS) { 
    Impl& surfImpl = SimTK_DYNAMIC_CAST_DEBUG<Impl&>(*impl);
    surfImpl.nearPointInS = stationInS; 
    return *this; 
}

// Points are provided in S frame but converted to coordinates.
CableObstacle::Surface& CableObstacle::Surface::
setContactPointHints(const Vec3& Phint_S, const Vec3& Qhint_S) { 
    Impl& surfImpl = SimTK_DYNAMIC_CAST_DEBUG<Impl&>(*impl);
    surfImpl.xPhint = Phint_S; // TODO: only works for implicit surfaces
    surfImpl.xQhint = Qhint_S; 
    return *this; 
}

//==============================================================================
//                   CABLE OBSTACLE :: SURFACE :: IMPL
//==============================================================================
void CableObstacle::Surface::Impl::initSurfacePath
   (const Vec3&        xP, 
    const Vec3&        xQ,
    const UnitVec3&    entryHint_S, // optional (use NaN)
    const UnitVec3&    exitHint_S,  // optional (use NaN)
    Geodesic&          path) const
{
    // TODO
    surface.initGeodesic(xP, xQ, nearPointInS, 
        GeodesicOptions(), path);
}

//------------------------------------------------------------------------------
//                       CALC SURFACE PATH ERROR
//------------------------------------------------------------------------------
Vec6 CableObstacle::Surface::Impl::calcSurfacePathError
   (const Geodesic& previous,
    const UnitVec3& eIn,    // from Qi-1 to Pi, in S
    const Vec3&     xP,     // coordinates of Pi
    const Vec3&     xQ,     // coordinates of Qi
    const UnitVec3& eOut,   // from Qi to Pi+1, in S
    Geodesic&       current) const
{
    surface.continueGeodesic(xP, xQ, previous, GeodesicOptions(), current);

    const UnitVec3& nP = current.getNormalP();
    const UnitVec3& nQ = current.getNormalQ();
    const UnitVec3& tP = current.getTangentP();
    const UnitVec3& tQ = current.getTangentQ();
    const UnitVec3& bP = current.getBinormalP();
    const UnitVec3& bQ = current.getBinormalQ();
    const Real      length = current.getLength();

    // Watch for backwards geodesic and flip tangent error conditions.
    Real signP = Real(~eIn *tP < 0 ? -1 : 1);
    Real signQ = Real(~eOut*tQ < 0 ? -1 : 1);
    //XXX 
    //signP = signQ = 1;

    // If length is very short, or geodesic is backwards, use path binormals
    // rather than geodesic binormals.
    const Real ShortLength = Real(1e-3);
    if (length <= ShortLength)
        cout << "==> Using short formulation for length=" << length << endl;
    if (signP < 0 || signQ < 0) { 
        cout << "==> Backwards geodesic with length=" << length 
             << " eIn*tP=" << ~eIn*tP << " eOut*tQ=" << ~eOut*tQ << endl;
        if (signP != signQ) {
            cout <<"    Mismatch eIn="<<eIn << " eOut="<<eOut<<endl;
            cout <<"             tP=" <<tP << " tQ="<<tQ<<endl;
        }
    }

    const Vec3 bbarP = length<=ShortLength ? eOut % nP : Vec3(bP);
    const Vec3 bbarQ = length<=ShortLength ? eIn  % nQ : Vec3(bQ);

    Vec6 err;
    err[0] = ~eIn*nP;   // tangent error in normal direction
    err[1] = ~eOut*nQ;
    
    err[2] = ~eIn*bbarP;   // tangent errors in geodesic direction
    err[3] = ~eOut*bbarQ;

    // These are the implicit surface errors forcing P and Q to lie on the
    // surface.
    err[4] = surface.calcSurfaceValue(xP);
    err[5] = surface.calcSurfaceValue(xQ);

    //cout << "  surfErr=" << err << endl;
    return err;
}



//------------------------------------------------------------------------------
//                CALC SURFACE PATH ERROR JACOBIAN ANALYTICALLY
//------------------------------------------------------------------------------
void CableObstacle::Surface::Impl::
calcSurfacePathErrorJacobianAnalytically
   (const UnitVec3& eIn,        // from Qi-1 to Pi, in S
    const Vec3&     xP,         // coordinates of Pi, in S
    const Vec3&     xQ,         // coordinates of Qi, in S
    const UnitVec3& eOut,       // from Qi to Pi+1, in S
    const Geodesic& current,    // Must have been computed from above params.
    Mat63&          DerrDentry, // 4x3 for parametric
    Mat63&          DerrDxP,    // 4x2       "
    Mat63&          DerrDxQ,    // 4x2       "
    Mat63&          DerrDexit)  // 4x3       "
    const
{
    const UnitVec3 nP = current.getNormalP();
    const UnitVec3 nQ = current.getNormalQ();
    const UnitVec3& tP = current.getTangentP();
    const UnitVec3& tQ = current.getTangentQ();
    const UnitVec3& bP = current.getBinormalP();
    const UnitVec3& bQ = current.getBinormalQ();
    // Watch for backwards geodesic and flip tangent error conditions.
    Real signP = Real(~eIn *tP < 0 ? -1 : 1);
    Real signQ = Real(~eOut*tQ < 0 ? -1 : 1);
    // XXX
    //signP = signQ = 1;

    const Real jP = current.getJacobiP(), 
               jQ = current.getJacobiQ();
    const Real jdP = current.getJacobiPDot(), 
               jdQ = current.getJacobiQDot();
    const Real tauP = current.getTorsionP(), 
               tauQ = current.getTorsionQ();
    const Real kappaP = current.getCurvatureP(), 
               kappaQ = current.getCurvatureQ();
    const Real muP = current.getBinormalCurvatureP(),
               muQ = current.getBinormalCurvatureQ();

    DerrDentry = Mat63( ~nP,    Row3(0), signP*~bP,   Row3(0), Row3(0), Row3(0) );
    DerrDexit  = Mat63( Row3(0), ~nQ,   Row3(0), signP*~bQ,    Row3(0), Row3(0) );
    const Vec3 gP = surface.calcSurfaceGradient(xP),
               gQ = surface.calcSurfaceGradient(xQ);
    const Mat33 HP = surface.calcSurfaceHessian(xP),
                HQ = surface.calcSurfaceHessian(xQ);
    const Real oojP = std::abs(jP) < SqrtEps ? Real(0) : 1/jP, 
               oojQ = std::abs(jQ) < SqrtEps ? Real(0) : 1/jQ;
    const Mat33 DnPDxP = (Mat33(1) - nP*~nP)*HP / (~gP*nP),
                DnQDxQ = (Mat33(1) - nQ*~nQ)*HQ / (~gQ*nQ);
    const Mat33 DbPDxP = -tauP*nP*~tP - (oojP*jdP*tP + muP*nP)*~bP,
                DbQDxQ = -tauQ*nQ*~tQ - (oojQ*jdQ*tQ + muQ*nQ)*~bQ;
    const Mat33 DbPDxQ = -oojQ*tP*~bQ,
                DbQDxP =  oojP*tQ*~bP;

    DerrDxP = Mat63( ~eIn  * DnPDxP,
                         Row3(0),
                     signP*(~eIn  * DbPDxP),
                     signQ*(~eOut * DbQDxP),
                          ~gP,
                         Row3(0)    );

    DerrDxQ = Mat63(     Row3(0),
                     ~eOut * DnQDxQ,
                     signP*(~eIn  * DbPDxQ),
                     signQ*(~eOut * DbQDxQ),
                         Row3(0),
                          ~gQ       );
}

//------------------------------------------------------------------------------
//                  CALC SURFACE KINEMATIC VELOCITY ERROR
//------------------------------------------------------------------------------
// Note: this must calculate the *negative* of the kinematic velocity term.
Vec6 CableObstacle::Surface::Impl::calcSurfaceNegKinematicVelocityError
   (const Geodesic& geodesic,
    const UnitVec3& eIn,        // from Qi-1 to Pi, in S
    const Vec3&     eInDot,     // time derivative of eIn, taken in S
    const Vec3&     xP,         // coordinates of Pi
    const Vec3&     xQ,         // coordinates of Qi
    const UnitVec3& eOut,       // from Qi to Pi+1, in S
    const Vec3&     eOutDot)    // time derivative of eOut, taken in S
    const
{
    Vec6 errdotK;

    const UnitVec3& nP = geodesic.getNormalP();
    const UnitVec3& nQ = geodesic.getNormalQ();
    const UnitVec3& tP = geodesic.getTangentP();
    const UnitVec3& tQ = geodesic.getTangentQ();
    const UnitVec3& bP = geodesic.getBinormalP();
    const UnitVec3& bQ = geodesic.getBinormalQ();

    // Watch for backwards geodesic and flip tangent error conditions.
    Real signP = Real(~eIn *tP < 0 ? -1 : 1);
    Real signQ = Real(~eOut*tQ < 0 ? -1 : 1);
    //XXX
    //signP = signQ = 1;

    // These must be the kinematic part of the time derivatives of the
    // error conditions that were reported for this obstacle during path
    // error calculation, but negated.
    errdotK[0] = - ~eInDot*nP; 
    errdotK[1] = - ~eOutDot*nQ;
    errdotK[2] = -signP*(~eInDot*bP);
    errdotK[3] = -signQ*(~eOutDot*bQ);
    // Implicit surface error is frozen since xP and xQ are.
    errdotK[4] = 0;
    errdotK[5] = 0;

    return errdotK;
}

// Numerical Jacobian calculation.

class SurfaceError : public Differentiator::JacobianFunction {
public:
    // Caution: we're holding references to surface and prevGeodesic.
    SurfaceError(const CableObstacle::Surface::Impl& surface,
                 const Geodesic& prevGeodesic,
                 Real accuracy)
    :   Differentiator::JacobianFunction(6, 12), surface(surface),
        prevGeodesic(prevGeodesic)
    {   setEstimatedAccuracy(accuracy); }

    int f(const Vector& x, Vector& fx) const override {
        UnitVec3 eIn, eOut;
        Vec3     xP, xQ;
        Geodesic newGeodesic; // throw away
        mapXToVecs(x, eIn, xP, xQ, eOut);
        Vec6 err = surface.calcSurfacePathError(prevGeodesic,
                             eIn, xP, xQ, eOut, newGeodesic);
        assert(fx.size() == 6);
        Vec6::updAs(&fx[0]) = err;
        return 0;
    }

    void mapXToVecs(const Vector& x, 
                    UnitVec3& eIn, Vec3& xP, Vec3& xQ, UnitVec3& eOut) const
    {
        assert(x.size()==12);
        //eIn = UnitVec3(Vec3::getAs(&x[0]), true); // no need to normalize
        eIn = UnitVec3(Vec3::getAs(&x[0]));
        xP = Vec3::getAs(&x[3]);
        xQ = Vec3::getAs(&x[6]);
        //eOut = UnitVec3(Vec3::getAs(&x[9]), true);
        eOut = UnitVec3(Vec3::getAs(&x[9]));
    }

    void mapVecsToX(const UnitVec3& eIn, const Vec3& xP, 
                    const Vec3& xQ, const UnitVec3& eOut,
                    Vector& x) const
    {
        assert(x.size()==12);
        Vec3::updAs(&x[0]) = eIn.asVec3();
        Vec3::updAs(&x[3]) = xP;
        Vec3::updAs(&x[6]) = xQ;
        Vec3::updAs(&x[9]) = eOut.asVec3();
    }
    
    // Map 6x12 matrix into four 6x3 blocks. We're assuming the Matrix is
    // compact and column-ordered.
    void mapMatrixToMats(const Matrix& J, 
        Mat63&          DerrDentry,
        Mat63&          DerrDxP,
        Mat63&          DerrDxQ,
        Mat63&          DerrDexit) const
    {
        DerrDentry = Mat63::getAs(&J(0,0));
        DerrDxP    = Mat63::getAs(&J(0,3));
        DerrDxQ    = Mat63::getAs(&J(0,6));
        DerrDexit  = Mat63::getAs(&J(0,9));
    }
private:
    const CableObstacle::Surface::Impl& surface;
    const Geodesic&                     prevGeodesic;
};

void CableObstacle::Surface::Impl::
calcSurfacePathErrorJacobianNumerically
   (const Geodesic& previous,
    const UnitVec3& eIn_S,
    const Vec3&     xP,
    const Vec3&     xQ,
    const UnitVec3& eOut_S,
    const Geodesic& current,    // Must have been computed from above params.
    Mat63&          DerrDentry, // 4x3 for parametric
    Mat63&          DerrDxP,    // 4x2       "
    Mat63&          DerrDxQ,    // 4x2       "
    Mat63&          DerrDexit)  // 4x3       "
    const
{
    SurfaceError jacFunction(*this, previous, Real(1e-8)); // TODO
    Differentiator diff(jacFunction);
    Vector x(12), err0(6);
    jacFunction.mapVecsToX(eIn_S,xP,xQ,eOut_S, x);
    jacFunction.f(x,err0);
    Matrix J(6,12);
    diff.calcJacobian(x, err0, J, Differentiator::CentralDifference);
    jacFunction.mapMatrixToMats(J, DerrDentry, DerrDxP, DerrDxQ, DerrDexit);
}
