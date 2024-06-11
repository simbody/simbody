/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2024 Authors.
 Authors: Pepijn van den Bos
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ----------------------------------------------------------------------------*/

#include "CableSpan_CurveSegment.h"
#include "CableSpan_Impl.h"

#include "simbody/internal/CableSpan.h"
#include "simbody/internal/MultibodySystem.h"

using namespace SimTK;

// Relevant aliases for convenience.
using InstanceEntry       = CurveSegment::InstanceEntry;
using LocalGeodesicSample = CurveSegment::LocalGeodesicSample;
using WrappingStatus      = CurveSegment::WrappingStatus;
using IntegratorTolerances = CurveSegment::IntegratorTolerances;
using ObstacleIndex        = CableSpan::ObstacleIndex;

//==============================================================================
//                      Cached Data Structures
//==============================================================================
// Begin anonymous namespace, limiting the data structures to this compilation
// unit.
namespace
{

//==============================================================================
//                         Struct MatrixWorkspace
//==============================================================================
/** This is a helper struct that is used by a CableSpan to compute the
Stage::Position level data, i.e. the spanned path.
Computing the spanned path envolves a Newton type iteration, for which several
matrices are computed.
After computing the path this data is no longer needed by the CableSpan. */
struct MatrixWorkspace
{
    // Number of degrees of freedom of a geodesic.
    static constexpr int GEODESIC_DOF = 4;

    // Number of path error constraints per curve segment.
    static constexpr int NUMBER_OF_CONSTRAINTS = 4;

    // Given the number of CurveSegments that are in contact with their
    // respective obstacle's surface, contruct a MatrixWorkspace of correct
    // dimensions.
    MatrixWorkspace(int nActive) : nCurves(nActive)
    {
        static constexpr int Q = GEODESIC_DOF;
        // 4 for the path error, and 1 for the weighting of the length.
        static constexpr int C = NUMBER_OF_CONSTRAINTS + 1;
        const int n            = nActive;

        lineSegments.resize(n + 1);
        pathErrorJacobian = Matrix(C * n, Q * n, 0.);
        pathCorrection    = Vector(Q * n, 0.);
        pathError         = Vector(C * n, 0.);
    }

    std::vector<LineSegment> lineSegments;

    Matrix pathErrorJacobian;
    Vector pathCorrection;
    Vector pathError;
    FactorQTZ inverse;
    int nCurves = -1;
};

//==============================================================================
//                      Struct PathSolverScratchData
//==============================================================================
// Cache entry for holding MatrixWorkspaces of different dimensions.
//
// CableSubsystem caches this data, such that all it's CableSpans can make use
// of it as a scratchpad when computing the Stage::Position level data.
class PathSolverScratchData
{
public:
    PathSolverScratchData() = default;

    // Get mutable access to a MatrixWorkspace of appropriate dimension.
    // Constructs requested MatrixWorkspace of requested dimension if not done previously.
    MatrixWorkspace& updOrInsert(int nActive)
    {
        SimTK_ASSERT1(nActive > 0,
            "PathSolverScratchData::updOrInsert()"
            "Number of obstacles in contact must be larger than zero (got %d)",
            nActive);

        // Construct all MatrixWorkspace's of requested dimension and lower.
        for (int i = matrixWorkspaces.size(); i < nActive; ++i) {
            matrixWorkspaces.emplace_back(i + 1);
        }

        // Return MatrixWorkspace of requested dimension.
        return matrixWorkspaces.at(nActive - 1);
    }

private:
    // MatrixWorkspaces of suitable dimensions for solving the cable path given
    // a number of obstacles in contact with the cable. The ith element in the
    // vector is used for solving a path with i+1 active obstacles.
    std::vector<MatrixWorkspace> matrixWorkspaces;
};

}

//==============================================================================
//                      Class CableSubsystem::Impl
//==============================================================================
/* CableSubsystem::Impl TODO desc. */
class CableSubsystem::Impl : public Subsystem::Guts
{
public:
    // TODO fix rule of five
    Impl() = default;

    void realizeTopology(State& state)
    {
        PathSolverScratchData cache{};
        m_CacheIx = allocateCacheEntry(
            state,
            Stage::Instance,
            Stage::Infinity,
            new Value<PathSolverScratchData>(cache));
    }

    CableSpanIndex adoptCable(CableSubsystem& subsystemHandle, CableSpan& cable)
    {
        invalidateSubsystemTopologyCache();

        CableSpanIndex cableIx(cables.size());
        cable.updImpl().setSubsystem(subsystemHandle, cableIx);
        cables.emplace_back(cable);
        return cableIx;
    }

    int getNumCables() const
    {
        return cables.size();
    }

    const CableSpan& getCable(CableSpanIndex index) const
    {
        SimTK_ERRCHK1_ALWAYS(
                &cables[index].getImpl().getSubsystem().getImpl() == this,
                "CableSubsystem::getCable", "Cable %d is not owned by this subsystem",
                index);
        SimTK_ERRCHK1_ALWAYS(
                cables[index].getImpl().getIndex() == index,
                "CableSubsystem::getCable",
                "Cable %d has an invalid index",
                index);
        return cables[index];
    }

    CableSpan& updCable(CableSpanIndex index)
    {
        SimTK_ERRCHK1_ALWAYS(
                &cables[index].getImpl().getSubsystem().getImpl() == this,
                "CableSubsystem::updCable", "Cable %d is not owned by this subsystem",
                index);
        SimTK_ERRCHK1_ALWAYS(
                cables[index].getImpl().getIndex() == index,
                "CableSubsystem::getCable",
                "Cable %d has an invalid index",
                index);
        return cables[index];
    }

    // Return the MultibodySystem which owns this WrappingPathSubsystem.
    const MultibodySystem& getMultibodySystem() const
    {
        return MultibodySystem::downcast(getSystem());
    }

    // TODO grab correct dimension here.
    PathSolverScratchData& updSolverData(const State& state) const
    {
        return Value<PathSolverScratchData>::updDowncast(updCacheEntry(state, m_CacheIx));
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    Impl* cloneImpl() const override
    {
        return new Impl(*this);
    }

    int realizeSubsystemTopologyImpl(State& state) const override
    {
        // Briefly allow writing into the Topology cache. After this the
        // Topology cache is const.
        Impl* mutableThis = const_cast<Impl*>(this);

        mutableThis->realizeTopology(state);
        for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
            CableSpan& path = mutableThis->updCable(ix);
            path.updImpl().realizeTopology(state);
        }

        return 0;
    }

    int calcDecorativeGeometryAndAppendImpl(
            const State& state,
            Stage stage,
            Array_<DecorativeGeometry>& decorations) const override
    {
        if (stage != Stage::Position) {
            return 0;
        }

        for (const CableSpan& cable : cables) {
            int returnValue = cable.getImpl().calcDecorativeGeometryAndAppend(
                    state,
                    stage,
                    decorations);
            if (returnValue != 0) {
                return returnValue;
            }
        }
        return 0;
    }

    // TOPOLOGY STATE
    Array_<CableSpan, CableSpanIndex> cables;

    CacheEntryIndex m_CacheIx;

    friend CableSubsystemTestHelper;
};

//==============================================================================
//                          CABLE SUBSYSTEM
//==============================================================================

bool CableSubsystem::isInstanceOf(const Subsystem& s)
{
    return Impl::isA(s.getSubsystemGuts());
}

const CableSubsystem& CableSubsystem::downcast(const Subsystem& s)
{
    assert(isInstanceOf(s));
    return static_cast<const CableSubsystem&>(s);
}

CableSubsystem& CableSubsystem::updDowncast(Subsystem& s)
{
    assert(isInstanceOf(s));
    return static_cast<CableSubsystem&>(s);
}

CableSubsystem::CableSubsystem()
{
    adoptSubsystemGuts(new Impl());
}

CableSubsystem::CableSubsystem(MultibodySystem& mbs)
{
    adoptSubsystemGuts(new Impl());
    mbs.adoptSubsystem(*this);
}

int CableSubsystem::getNumCables() const
{
    return getImpl().getNumCables();
}

const CableSpan& CableSubsystem::getCable(CableSpanIndex ix) const
{
    return getImpl().getCable(ix);
}

CableSpan& CableSubsystem::updCable(CableSpanIndex ix)
{
    return updImpl().updCable(ix);
}

const CableSubsystem::Impl& CableSubsystem::getImpl() const
{
    return SimTK_DYNAMIC_CAST_DEBUG<const Impl&>(getSubsystemGuts());
}

CableSubsystem::Impl& CableSubsystem::updImpl()
{
    return SimTK_DYNAMIC_CAST_DEBUG<Impl&>(updSubsystemGuts());
}

//==============================================================================
//                            Frenet Frame Helpers
//==============================================================================

namespace
{

// Define the frenet frame as the position and orientation with:
// - position on the geodesic,
// - tangent along X direction,
// - surface normal along Y direction,
// - binormal along Z direction.
using FrenetFrame                        = Transform;
static const CoordinateAxis TangentAxis  = CoordinateAxis::XCoordinateAxis();
static const CoordinateAxis NormalAxis   = CoordinateAxis::YCoordinateAxis();
static const CoordinateAxis BinormalAxis = CoordinateAxis::ZCoordinateAxis();

const UnitVec3& getAxis(const FrenetFrame& X, CoordinateAxis axis)
{
    return X.R().getAxisUnitVec(axis);
}

const UnitVec3& getTangent(const FrenetFrame& X)
{
    return getAxis(X, TangentAxis);
}

const UnitVec3& getNormal(const FrenetFrame& X)
{
    return getAxis(X, NormalAxis);
}

const UnitVec3& getBinormal(const FrenetFrame& X)
{
    return getAxis(X, BinormalAxis);
}

} // namespace

//==============================================================================
//              Contact Geometry Helpers For Sign Inversion
//==============================================================================

// This section contains helper functions for flipping the sign of the output
// obtained from ContactGeometry.
namespace
{

// Flip sign of curvature, such that tangentDot = curvature * normal
Real calcSurfaceCurvature(
    const ContactGeometry& geometry,
    const FrenetFrame& X_S,
    CoordinateAxis direction)
{
    return -geometry.calcSurfaceCurvatureInDirection(
        X_S.p(),
        getAxis(X_S, direction));
}

// Helper function for flipping the sign such that positive value is outside of
// the surface.
Real calcSurfaceValue(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceValue(point_S);
}
// Flipping sign...
Vec3 calcSurfaceGradient(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceGradient(point_S);
}
// Flipping sign...
Mat33 calcSurfaceHessian(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceHessian(point_S);
}

} // namespace

//==============================================================================
//                                CABLE SPAN
//==============================================================================

// Allocate a new default-constructed CableSpan::Impl.
CableSpan::CableSpan() : m_Impl(new Impl()) {}

CableSpan::CableSpan(
    CableSubsystem& subsystem,
    MobilizedBodyIndex originBody,
    const Vec3& defaultOriginPoint,
    MobilizedBodyIndex terminationBody,
    const Vec3& defaultTerminationPoint) :
    m_Impl(std::shared_ptr<Impl>(new Impl(
        originBody,
        defaultOriginPoint,
        terminationBody,
        defaultTerminationPoint)))
{
    subsystem.updImpl().adoptCable(subsystem, *this);
}

//------------------------------------------------------------------------------

ObstacleIndex CableSpan::addSurfaceObstacle(
    MobilizedBodyIndex body,
    const Transform& X_BS,
    const ContactGeometry& geometry,
    Vec3 contactPointHint)
{
    return updImpl().addSurfaceObstacle(body, X_BS, geometry, contactPointHint);
}

int CableSpan::getNumSurfaceObstacles() const
{
    return getImpl().getNumCurveSegments();
}

const MobilizedBody& CableSpan::getObstacleMobilizedBody(ObstacleIndex ix) const
{
    return getImpl().getCurveSegment(ix).getMobilizedBody();
}

const MobilizedBodyIndex& CableSpan::getObstacleMobilizedBodyIndex(
    ObstacleIndex ix) const
{
    return getImpl().getCurveSegment(ix).getMobilizedBodyIndex();
}

void CableSpan::setObstacleMobilizedBodyIndex(
    ObstacleIndex ix,
    MobilizedBodyIndex body)
{
    updImpl().updCurveSegment(ix).setMobilizedBodyIndex(body);
}

const Transform& CableSpan::getObstacleXformSurfaceToBody(
    ObstacleIndex ix) const
{
    return getImpl().getCurveSegment(ix).getXformSurfaceToBody();
}
void CableSpan::setObstacleXformSurfaceToBody(ObstacleIndex ix, const Transform& X_BS)
{
    updImpl().updCurveSegment(ix).setXformSurfaceToBody(std::move(X_BS));
}

const ContactGeometry& CableSpan::getObstacleContactGeometry(
    ObstacleIndex ix) const
{
    return getImpl().getCurveSegment(ix).getContactGeometry();
}
void CableSpan::setObstacleContactGeometry(
    ObstacleIndex ix,
    ContactGeometry geometry)
{
    updImpl().updCurveSegment(ix).setContactGeometry(geometry);
}

Vec3 CableSpan::getObstacleInitialContactPointHint(ObstacleIndex ix) const
{
    return getImpl().getCurveSegment(ix).getContactPointHint();
}

void CableSpan::setObstacleInitialContactPointHint(
    ObstacleIndex ix,
    Vec3 initialContactPointHint)
{
    updImpl().updCurveSegment(ix).setContactPointHint(initialContactPointHint);
}

//------------------------------------------------------------------------------

WrappingStatus CableSpan::getObstacleWrappingStatus(
    const State& state,
    ObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).getInstanceEntry(state).status;
}

Real CableSpan::getCurveSegmentLength(const State& state, ObstacleIndex ix)
    const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).getInstanceEntry(state).length;
}

int CableSpan::calcCurveSegmentPathPoints(
    const State& state,
    ObstacleIndex ix,
    int nPoints,
    const std::function<void(Real length, Vec3 point)>& sink) const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).calcPathPoints(state, sink, nPoints);
}

// This is useful for debugging and visualization.
int CableSpan::calcCurveSegmentPathPointsAndTangents(
    const State& state,
    ObstacleIndex ix,
    int nPoints,
    const std::function<void(Real length, Vec3 point, UnitVec3 tangent)>& sink) const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).calcPathPointsAndTangents(
        state,
        sink,
        nPoints);
}

int CableSpan::getCurveSegmentNumberOfIntegratorStepsTaken(
    const State& state,
    ObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).getInstanceEntry(state).samples.size();
}

Real CableSpan::getCurveSegmentInitialIntegratorStepSize(
    const State& state,
    ObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl()
        .getCurveSegment(ix)
        .getInstanceEntry(state)
        .integratorInitialStepSize;
}

//------------------------------------------------------------------------------

Real CableSpan::getSurfaceConstraintTolerance() const
{
    return getImpl().getSurfaceConstraintTolerance();
}

void CableSpan::setSurfaceConstraintTolerance(Real accuracy)
{
    updImpl().setSurfaceConstraintTolerance(accuracy);
}

int CableSpan::getSurfaceProjectionMaxIterations() const
{
    return getImpl().getSurfaceProjectionMaxIter();
}

void CableSpan::setSurfaceProjectionMaxIterations(int maxIterations)
{
    updImpl().setSurfaceProjectionMaxIter(maxIterations);
}

Real CableSpan::getIntegratorAccuracy() const
{
    return getImpl().getIntegratorAccuracy();
}

void CableSpan::setIntegratorAccuracy(Real accuracy)
{
    updImpl().setIntegratorAccuracy(accuracy);
}

int CableSpan::getSolverMaxIter() const
{
    return getImpl().getSolverMaxIterations();
}

void CableSpan::setSolverMaxIter(int maxIterations)
{
    updImpl().setSolverMaxIterations(maxIterations);
}

Real CableSpan::getPathErrorAccuracy() const
{
    return getImpl().getPathAccuracy();
}

void CableSpan::setPathErrorAccuracy(Real accuracy)
{
    updImpl().setPathAccuracy(accuracy);
}

Real CableSpan::getMaxRadialStepInDegrees() const
{
    return getImpl().getMaxRadialStepInDegrees();
}

void CableSpan::setMaxRadialStepInDegrees(Real maxStepInDeg)
{
    updImpl().setMaxRadialStepInDegrees(maxStepInDeg);
}

//------------------------------------------------------------------------------

Real CableSpan::getLength(const State& s) const
{
    return getImpl().getPosInfo(s).cableLength;
}

Real CableSpan::getLengthDot(const State& s) const
{
    return getImpl().getVelInfo(s).lengthDot;
}

void CableSpan::applyBodyForces(
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    return getImpl().applyBodyForces(s, tension, bodyForcesInG);
}

int CableSpan::calcPathPoints(
    const State& state,
    Real lengthIncrement,
    const std::function<void(Real length, Vec3 point_G)>& sink) const
{
    return getImpl().calcPathPoints(state, lengthIncrement, sink);
}

Real CableSpan::calcCablePower(const State& state, Real tension) const
{
    return getImpl().calcCablePower(state, tension);
}

int CableSpan::getNumSolverIter(const State& state) const
{
    return getImpl().getPosInfo(state).loopIter;
}

Real CableSpan::getMaxPathError(const State& state) const
{
    return getImpl().getPosInfo(state).pathError;
}

void CableSpan::storeCurrentPath(State& state) const
{
    for (ObstacleIndex ix(0); ix < getImpl().getNumCurveSegments(); ++ix)
    {
        getImpl().getCurveSegment(ix).storeCurrentPath(state);
    }
}

//------------------------------------------------------------------------------
//                              CABLE SPAN IMPL
//------------------------------------------------------------------------------

void CableSpan::Impl::realizeTopology(State& s)
{
    for (CurveSegment& segment : m_CurveSegments) {
        segment.realizeTopology(s);
    }

    PosInfo ppe{};
    m_PosInfoIx = updSubsystem().allocateCacheEntry(
        s,
        Stage::Position,
        Stage::Infinity,
        new Value<PosInfo>(ppe));

    getSubsystem().markCacheValueNotRealized(s, m_PosInfoIx);

    VelInfo velInfo{};
    m_VelInfoIx = updSubsystem().allocateCacheEntry(
        s,
        Stage::Velocity,
        Stage::Infinity,
        new Value<VelInfo>(velInfo));

    getSubsystem().markCacheValueNotRealized(s, m_VelInfoIx);
}

void CableSpan::Impl::realizePosition(const State& s) const
{
    if (getSubsystem().isCacheValueRealized(s, m_PosInfoIx)) {
        return;
    }
    calcPosInfo(s, updPosInfo(s));
    getSubsystem().markCacheValueRealized(s, m_PosInfoIx);
}

void CableSpan::Impl::realizeVelocity(const State& s) const
{
    realizePosition(s);
    if (getSubsystem().isCacheValueRealized(s, m_VelInfoIx)) {
        return;
    }
    calcVelInfo(s, updVelInfo(s));
    getSubsystem().markCacheValueRealized(s, m_VelInfoIx);
}

void CableSpan::Impl::invalidatePositionLevelCache(const State& s) const
{
    for (const CurveSegment& curve : m_CurveSegments) {
        curve.invalidatePosEntry(s);
    }
    getSubsystem().markCacheValueNotRealized(s, m_PosInfoIx);
}

const CableSpan::Impl::PosInfo& CableSpan::Impl::getPosInfo(
    const State& s) const
{
    realizePosition(s);
    return Value<PosInfo>::downcast(
        getSubsystem().getCacheEntry(s, m_PosInfoIx));
}

CableSpan::Impl::PosInfo& CableSpan::Impl::updPosInfo(const State& s) const
{
    return Value<PosInfo>::updDowncast(
        getSubsystem().updCacheEntry(s, m_PosInfoIx));
}

const CableSpan::Impl::VelInfo& CableSpan::Impl::getVelInfo(
    const State& s) const
{
    realizeVelocity(s);
    return Value<VelInfo>::downcast(
        getSubsystem().getCacheEntry(s, m_VelInfoIx));
}

CableSpan::Impl::VelInfo& CableSpan::Impl::updVelInfo(const State& s) const
{
    return Value<VelInfo>::updDowncast(
        getSubsystem().updCacheEntry(s, m_VelInfoIx));
}

//------------------------------------------------------------------------------

const Mobod& CableSpan::Impl::getOriginBody() const
{
    return getSubsystem()
        .getImpl()
        .getMultibodySystem()
        .getMatterSubsystem()
        .getMobilizedBody(m_OriginBody);
}

const Mobod& CableSpan::Impl::getTerminationBody() const
{
    return getSubsystem()
        .getImpl()
        .getMultibodySystem()
        .getMatterSubsystem()
        .getMobilizedBody(m_TerminationBody);
}

//------------------------------------------------------------------------------

int CableSpan::Impl::countActive(const State& s) const
{
    int count = 0;
    for (const CurveSegment& segment : m_CurveSegments) {
        if (segment.isInContactWithSurface(s)) {
            ++count;
        }
    }
    return count;
}

int CableSpan::Impl::calcPathPoints(
    const State& state,
    Real maxLengthIncrement,
    const std::function<void(Real length, Vec3 point_G)>& sink) const
{
    // Count number of points written.
    int count = 0;

    // Write the initial point.
    const PosInfo& pos = getPosInfo(state);
    sink(0., pos.originPoint_G);
    ++count;

    // Write points along each of the curves.
    for (const CurveSegment& curve : m_CurveSegments) {

        Real curveLength = curve.getInstanceEntry(state).length;
        int nPoints =
            static_cast<int>(std::abs(curveLength / maxLengthIncrement) + 1.);
        nPoints = nPoints < 2 ? 2 : nPoints;

        count += curve.calcPathPoints(state, sink, nPoints);
    }

    // Write the termination point.
    sink(pos.cableLength, pos.terminationPoint_G);
    ++count;

    // Return number of points written.
    return count;
}

//------------------------------------------------------------------------------

const CurveSegment* CableSpan::Impl::findPrevActiveCurveSegment(
    const State& s,
    ObstacleIndex ix) const
{
    for (int i = ix - 1; i >= 0; --i) {
        // Find the active segment before the current.
        if (m_CurveSegments.at(ObstacleIndex(i))
                .isInContactWithSurface(s)) {
            return &m_CurveSegments.at(ObstacleIndex(i));
        }
    }
    return nullptr;
}

const CurveSegment* CableSpan::Impl::findNextActiveCurveSegment(
    const State& s,
    ObstacleIndex ix) const
{
    // Find the active segment after the current.
    for (int i = ix + 1; i < m_CurveSegments.size(); ++i) {
        if (m_CurveSegments.at(ObstacleIndex(i))
                .isInContactWithSurface(s)) {
            return &m_CurveSegments.at(ObstacleIndex(i));
        }
    }
    return nullptr;
}

Vec3 CableSpan::Impl::findPrevPoint(const State& s, ObstacleIndex ix) const
{
    // Check if there is a curve segment preceding given obstacle.
    const CurveSegment* prevCurve = findPrevActiveCurveSegment(s, ix);
    if (prevCurve) {
        // If so, the previous point is the final contact point of the previous
        // curve.
        const Vec3& finalContactPoint_S =
            prevCurve->getInstanceEntry(s).X_SQ.p();
        // Transform the contact point to ground.
        const Transform& X_GS = prevCurve->calcSurfaceFrameInGround(s);
        return X_GS.shiftFrameStationToBase(finalContactPoint_S);
    }
    // There are no curve segments before given obstacle: the previous point is
    // the path's origin point.
    return getOriginBody().getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
}

Vec3 CableSpan::Impl::findNextPoint(const State& s, ObstacleIndex ix) const
{
    // Check if there is a curve segment after given obstacle.
    const CurveSegment* nextCurve = findNextActiveCurveSegment(s, ix);
    if (nextCurve) {
        // If so, the next point is the initial contact point of the next
        // curve.
        const Vec3& initialContactPoint_S =
            nextCurve->getInstanceEntry(s).X_SP.p();
        // Transform the contact point to ground.
        const Transform& X_GS = nextCurve->calcSurfaceFrameInGround(s);
        return X_GS.shiftFrameStationToBase(initialContactPoint_S);
    };
    // There are no curve segments following given obstacle: the next point is
    // the path's termination point.
    return getTerminationBody()
        .getBodyTransform(s)
        .shiftFrameStationToBase(m_TerminationPoint);
}

//==============================================================================
//                      CABLE FORCE COMPUTATIONS
//==============================================================================

namespace
{

void calcUnitForceExertedByCurve(
    const CurveSegment& curve,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CurveSegment::PosEntry& ppe = curve.getPosInfo(s);
    const Vec3& x_GB = curve.getMobilizedBody().getBodyOriginLocation(s);

    // Contact point moment arms in ground.
    const Vec3 r_P = ppe.X_GP.p() - x_GB;
    const Vec3 r_Q = ppe.X_GQ.p() - x_GB;

    // Tangent directions at contact points in ground.
    const UnitVec3& t_P = getTangent(ppe.X_GP);
    const UnitVec3& t_Q = getTangent(ppe.X_GQ);

    unitForce_G[0] = r_Q % t_Q - r_P % t_P;
    unitForce_G[1] = t_Q - t_P;
}

void calcUnitForceAtCableOrigin(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    // Origin contact point moment arm in ground.
    const CableSpan::Impl::PosInfo ppe = cable.getPosInfo(s);
    const Vec3& arm_G =
        ppe.originPoint_G - cable.getOriginBody().getBodyOriginLocation(s);

    unitForce_G[0] = -arm_G % ppe.originTangent_G;
    unitForce_G[1] = -Vec3(ppe.originTangent_G);
}

void calcUnitForceAtCableTermination(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CableSpan::Impl::PosInfo ppe = cable.getPosInfo(s);
    const Vec3& arm_G                  = ppe.terminationPoint_G -
                        cable.getTerminationBody().getBodyOriginLocation(s);

    unitForce_G[0] = -arm_G % ppe.terminationTangent_G;
    unitForce_G[1] = -Vec3(ppe.terminationTangent_G);
}

} // namespace

void CableSpan::Impl::applyBodyForces(
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    realizePosition(s);

    if (tension < 0.) {
        return;
    }

    SpatialVec unitForce_G;

    // Force applied at cable origin point.
    {
        calcUnitForceAtCableOrigin(*this, s, unitForce_G);
        getOriginBody().applyBodyForce(s, unitForce_G * tension, bodyForcesInG);
    }

    // Forces applied to each obstacle body.
    for (const CurveSegment& curve : m_CurveSegments) {
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }

        calcUnitForceExertedByCurve(curve, s, unitForce_G);
        curve.getMobilizedBody().applyBodyForce(
            s,
            unitForce_G * tension,
            bodyForcesInG);
    }

    // Force applied at cable termination point.
    {
        calcUnitForceAtCableTermination(*this, s, unitForce_G);
        getTerminationBody().applyBodyForce(
            s,
            unitForce_G * tension,
            bodyForcesInG);
    }
}

Real CableSpan::Impl::calcCablePower(const State& state, Real tension) const
{
    realizePosition(state);

    if (tension < 0.) {
        return 0.;
    }

    SpatialVec unitForce_G;

    Real unitPower = 0.;

    {
        calcUnitForceAtCableOrigin(*this, state, unitForce_G);
        SpatialVec v = getOriginBody().getBodyVelocity(state);
        unitPower += ~unitForce_G * v;
    }

    for (const CurveSegment& curve : m_CurveSegments) {
        if (!curve.isInContactWithSurface(state)) {
            continue;
        }

        calcUnitForceExertedByCurve(curve, state, unitForce_G);
        SpatialVec v = curve.getMobilizedBody().getBodyVelocity(state);
        unitPower += ~unitForce_G * v;
    }

    {
        calcUnitForceAtCableTermination(*this, state, unitForce_G);
        SpatialVec v = getTerminationBody().getBodyVelocity(state);
        unitPower += ~unitForce_G * v;
    }

    return unitPower * tension;
}

//==============================================================================
//                             CABLE VISUALIZATION
//==============================================================================

namespace
{

// Helper for drawing a pretty curve segment.
void calcCurveDecorativeGeometryAndAppend(
    const CurveSegment& curve,
    const State& s,
    Array_<DecorativeGeometry>& decorations)
{
    constexpr Real c_InactiveOpacity = 0.25;

    const bool isInContactWithSurface = curve.isInContactWithSurface(s);
    const CurveSegment::PosEntry& ppe = curve.getPosInfo(s);

    // Draw the surface (TODO since we do not own it, should it be done here at
    // all?).
    {
        DecorativeGeometry geo = curve.getDecoration(); // TODO clone it?
        const Transform& X_GS  = ppe.X_GS;
        const Transform& X_SD  = geo.getTransform();
        // Inactive surfaces are dimmed.
        const Vec3 color = isInContactWithSurface
                               ? geo.getColor()
                               : geo.getColor() * c_InactiveOpacity;

        decorations.push_back(geo.setTransform(X_GS * X_SD).setColor(color));
    }

    // Check wrapping status to see if there is a curve to draw.
    if (!isInContactWithSurface) {
        return;
    }

    // Draw the curve segment as straight lines between prevPoint and nextPoint.
    {
        // Helper function for drawing line between prevPoint and nextPoint.
        Vec3 prevPoint{NaN};
        std::function<void(Real, Vec3)> drawLine = [&](Real length,
                                                       Vec3 nextPoint) {
            if (length != 0.) {
                decorations.push_back(DecorativeLine(prevPoint, nextPoint)
                                          .setColor(Purple)
                                          .setLineThickness(3));
            }
            prevPoint = nextPoint;
        };

        // TODO for debugging: Choose to show the original intergrator's path
        // points or the resampled path points.
        constexpr bool doResampling = false;
        if (doResampling) {
            curve.calcPathPoints(s, drawLine);
        } else {
            // Use points from GeodesicIntegrator directly.
            const std::vector<LocalGeodesicSample>& samples =
                curve.getInstanceEntry(s).samples;
            for (const LocalGeodesicSample& sample : samples) {
                drawLine(
                    sample.length,
                    ppe.X_GS.shiftFrameStationToBase(sample.frame.p()));
            }
        }
    }

    // Draw the Frenet frame at curve start and end.
    // TODO this is for debugging should be removed.
    {
        constexpr int FRENET_FRAME_LINE_THICKNESS = 5;
        constexpr Real FRENET_FRAME_LINE_LENGTH   = 0.5;

        const std::array<CoordinateAxis, 3> axes = {
            TangentAxis,
            NormalAxis,
            BinormalAxis};
        const std::array<Vec3, 3> colors = {Red, Green, Blue};

        std::function<void(const FrenetFrame& K)> DrawFrenetFrame =
            [&](const FrenetFrame& K) {
                for (int i = 0; i < 3; ++i) {
                    decorations.push_back(
                        DecorativeLine(
                            K.p(),
                            K.p() + FRENET_FRAME_LINE_LENGTH *
                                        K.R().getAxisUnitVec(axes.at(i)))
                            .setColor(colors.at(i))
                            .setLineThickness(FRENET_FRAME_LINE_THICKNESS));
                }
            };

        DrawFrenetFrame(ppe.X_GP);
        DrawFrenetFrame(ppe.X_GQ);
    }

    // Draw the initial contact point hint using a yellow line.
    // TODO this is for debugging should be removed.
    {
        constexpr int c_HintLineThickness = 2;

        const Transform& X_GS = ppe.X_GS;

        decorations.push_back(
            DecorativeLine(
                ppe.X_GS.p(),
                ppe.X_GS.shiftFrameStationToBase(curve.getContactPointHint()))
                .setColor(Yellow)
                .setLineThickness(c_HintLineThickness));
    }
}

} // namespace

int CableSpan::Impl::calcDecorativeGeometryAndAppend(
    const State& s,
    Stage stage,
    Array_<DecorativeGeometry>& decorations) const
{
    constexpr int c_LineThickness = 3;

    const PosInfo& ppe = getPosInfo(s);
    Vec3 prevPoint     = ppe.originPoint_G;

    for (const CurveSegment& curve : m_CurveSegments) {
        if (curve.isInContactWithSurface(s)) {
            const CurveSegment::PosEntry cppe = curve.getPosInfo(s);

            const Vec3 nextPoint = cppe.X_GP.p();
            decorations.push_back(DecorativeLine(prevPoint, nextPoint)
                                      .setColor(Purple)
                                      .setLineThickness(c_LineThickness));
            prevPoint = cppe.X_GQ.p();
        }

        calcCurveDecorativeGeometryAndAppend(curve, s, decorations);
    }

    // TODO choose colors.
    decorations.push_back(DecorativeLine(prevPoint, ppe.terminationPoint_G)
                              .setColor(Purple)
                              .setLineThickness(3));

    return 0;
}

//==============================================================================
//                             PATH ERROR VECTOR
//==============================================================================

namespace
{

// Helper for computing a single element of the path error vector.
Real calcPathError(const LineSegment& e, const Rotation& R, CoordinateAxis axis)
{
    return dot(e.direction, R.getAxisUnitVec(axis));
}

// Iterate over the cable's curve segments, and call provided function for
// those that are in contact with the obstacle's surfacce.
int forEachActiveCurveSegment(
    const CableSpan::Impl& cable,
    const State& s,
    const std::function<void(const CurveSegment& curve)>& callMe)
{
    int nActive = 0;
    for (ObstacleIndex ix(0); ix < cable.getNumCurveSegments(); ++ix) {
        const CurveSegment& curve = cable.getCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }
        callMe(curve);
        ++nActive;
    }
    return nActive;
}

// Compute the straight line segments in between the obstacles.
void calcLineSegments(
    const CableSpan::Impl& cable,
    const State& s,
    Vec3 pathOriginPoint,
    Vec3 pathTerminationPoint,
    std::vector<LineSegment>& lines)
{
    lines.clear();
    Vec3 prevPathPoint(pathOriginPoint);

    const int nActive =
        forEachActiveCurveSegment(cable, s, [&](const CurveSegment& curve) {
            const CurveSegment::PosEntry& curvePE = curve.getPosInfo(s);

            // Compute the line segment from the previous point to the
            // first curve contact point.
            lines.emplace_back(prevPathPoint, curvePE.getFirstContactPoint());

            // The next line segment will start at the curve's final contact
            // point.
            prevPathPoint = curvePE.getFinalContactPoint();
        });

    // Compute the last line segment.
    lines.emplace_back(prevPathPoint, pathTerminationPoint);

    // Sanity check.
    if (lines.size() != nActive + 1) {
        throw std::runtime_error(
            "Number of line segments does not match number of curve segments");
    }
}

Real calcTotalCableLength(
    const CableSpan::Impl& cable,
    const State& s,
    const std::vector<LineSegment>& lines)
{
    Real totalCableLength = 0.;

    // Add length of all line segments.
    for (const LineSegment& line : lines) {
        totalCableLength += line.length;
    }

    // Add length of each curve segment.
    const int nActive =
        forEachActiveCurveSegment(cable, s, [&](const CurveSegment& curve) {
            // Update the total cable length.
            totalCableLength += curve.getInstanceEntry(s).length;
        });

    // Sanity check.
    if (lines.size() != nActive + 1) {
        throw std::runtime_error(
            "Number of line segments does not match number of curve segments");
    }

    return totalCableLength;
}

} // namespace

template <size_t N>
void CableSpan::Impl::calcPathErrorVector(
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError) const
{
    int lineIx = 0;
    int row = -1;

    // Reset path error vector to zero.
    pathError *= 0;

    forEachActiveCurveSegment(
        *this,
        s,
        // Compute the path error for each curved segment.
        [&](const CurveSegment& curve) {
            const CurveSegment::PosEntry& ppe = curve.getPosInfo(s);
            // Compute path error at first contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++row) =
                    calcPathError(lines.at(lineIx), ppe.X_GP.R(), axis);
            }
            ++lineIx;
            // Compute path error at final contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++row) =
                    calcPathError(lines.at(lineIx), ppe.X_GQ.R(), axis);
            }
        });
}

//==============================================================================
//                             PATH ERROR JACOBIAN
//==============================================================================

// This section contains helper functions for computing the contribution of each
// curve segment to the total path error jacobian.
namespace
{

Vec4 calcJacobianOfPrevPathError(
    const CurveSegment& curve,
    const State& state,
    const LineSegment& line,
    const UnitVec3& axis)
{
    const Transform& X_GP = curve.getPosInfo(state).X_GP;

    // Compute partial derivative of path error to first contact point position
    // (x_PS), represented in the Frenet frame.
    Vec3 dErrDx = axis - line.direction * dot(line.direction, axis);
    dErrDx      = X_GP.RInv() * (dErrDx / line.length);

    return {dErrDx[0], dErrDx[2], 0., 0.};
}

Vec4 calcJacobianOfNextPathError(
    const CurveSegment& curve,
    const State& s,
    const LineSegment& line,
    const UnitVec3& axis)
{
    const Transform& X_GQ = curve.getPosInfo(s).X_GQ;

    const InstanceEntry& ie = curve.getInstanceEntry(s);
    const Real a            = ie.jacobi_Q[0];
    const Real r            = ie.jacobi_Q[1];

    // Partial derivative of path error to contact point position (x_QS),
    // represented in the frenet frame.
    Vec3 dErrDx = axis - line.direction * dot(line.direction, axis);
    dErrDx      = X_GQ.RInv() * (-dErrDx / line.length);

    return {dErrDx[0], dErrDx[2] * a, dErrDx[2] * r, dErrDx[0]};
}

Vec4 calcJacobianOfPathErrorAtP(
    const CurveSegment& curve,
    const State& s,
    const LineSegment& line,
    const UnitVec3& axis)
{
    const Transform& X_GP = curve.getPosInfo(s).X_GP;

    const InstanceEntry& ie = curve.getInstanceEntry(s);

    const Real tau = ie.torsion_P;
    const Real kt  = ie.curvatures_P[0];
    const Real kb  = ie.curvatures_P[1];

    // Compute partial derivative of path error to the natural geodesic
    // corrections. For derivation see TODO.
    const Vec3 y = X_GP.RInv() * cross(axis, line.direction);
    return Vec4{tau * y[0] + kt * y[2], -kb * y[0] - tau * y[2], -y[1], 0.} +
           calcJacobianOfPrevPathError(curve, s, line, axis);
}

Vec4 calcJacobianOfPathErrorAtQ(
    const CurveSegment& curve,
    const State& s,
    const LineSegment& line,
    const UnitVec3& axis)
{
    const Transform& X_GQ = curve.getPosInfo(s).X_GQ;

    const InstanceEntry& ie = curve.getInstanceEntry(s);

    const Real tau = ie.torsion_Q;
    const Real kt  = ie.curvatures_Q[0];
    const Real kb  = ie.curvatures_Q[1];

    const Real a = ie.jacobi_Q[0];
    const Real r = ie.jacobi_Q[1];

    const Real aDot = ie.jacobiDot_Q[0];
    const Real rDot = ie.jacobiDot_Q[1];

    // Partial derivative of path error to the given axis, represented in the
    // frenet frame.
    const Vec3 y = X_GQ.RInv() * cross(axis, line.direction);

    return Vec4{
               tau * y[0] + kt * y[2],
               -a * kb * y[0] - aDot * y[1] - a * tau * y[2],
               -r * kb * y[0] - rDot * y[1] - r * tau * y[2],
               tau * y[0] + kt * y[2]} +
           calcJacobianOfNextPathError(curve, s, line, axis);
}

} // namespace

template <size_t N>
void CableSpan::Impl::calcPathErrorJacobian(
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J) const
{
    constexpr int Nq = MatrixWorkspace::GEODESIC_DOF;

    // TODO perhaps just not make method static.
    const int n = lines.size() - 1;
    J *= 0.;

    SimTK_ASSERT(
        J.cols() == n * Nq,
        "Invalid number of columns in jacobian matrix");

    int row = 0;
    int col = 0;

    auto AddBlock = [&](const Vec4& block, int colOffset = 0) {
        for (int ix = 0; ix < 4; ++ix) {
            J.updElt(row, col + colOffset + ix) += block[ix];
        }
    };

    auto linesIt = lines.begin();
    for (ObstacleIndex ix(0); ix < getNumCurveSegments(); ++ix) {
        const CurveSegment& curve = getCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }

        const LineSegment& l_P = *linesIt;
        const LineSegment& l_Q = *++linesIt;

        const CurveSegment* prev = findPrevActiveCurveSegment(s, ix);
        const CurveSegment* next = findNextActiveCurveSegment(s, ix);

        const CurveSegment::PosEntry& ppe = curve.getPosInfo(s);
        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_P = getAxis(ppe.X_GP, axis);

            AddBlock(calcJacobianOfPathErrorAtP(curve, s, l_P, a_P));

            if (prev) {
                AddBlock(calcJacobianOfNextPathError(*prev, s, l_P, a_P), -Nq);
            }
            ++row;
        }

        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_Q = getAxis(ppe.X_GQ, axis);

            AddBlock(calcJacobianOfPathErrorAtQ(curve, s, l_Q, a_Q));

            if (next) {
                AddBlock(calcJacobianOfPrevPathError(*next, s, l_Q, a_Q), Nq);
            }
            ++row;
        }

        col += Nq;
    };
}

//==============================================================================
//                      SOLVING FOR GEODESIC CORRECTIONS
//==============================================================================

namespace
{

constexpr int GEODESIC_DOF          = MatrixWorkspace::GEODESIC_DOF;
constexpr int NUMBER_OF_CONSTRAINTS = MatrixWorkspace::NUMBER_OF_CONSTRAINTS;

// Natural geodesic correction vector.
//
// The elements are (in order):
// - Tangential correction of initial contact point,
// - Binormal correction of initial contact point,
// - Directional correction of initial contact tangent,
// - Lengthening correction of geodesic
using Correction = Vec<GEODESIC_DOF>;

// Solve for the geodesic corrections by attempting to set the path error to
// zero. We call this after having filled in the pathError vector and pathError
// jacobian in the MatrixWorkspace. The result is a vector of Corrections for each
// curve.
void calcPathCorrections(MatrixWorkspace& data, Real weight)
{
    // TODO add explanation...
    // Add a cost to changing the length.
    for (int i = 0; i < data.nCurves; ++i) {
        int r = data.nCurves * NUMBER_OF_CONSTRAINTS + i;
        int c = GEODESIC_DOF * (i + 1) - 1;
        data.pathErrorJacobian.set(r, c, weight);
    }

    data.inverse = data.pathErrorJacobian;
    data.inverse.solve(data.pathError, data.pathCorrection);
    data.pathCorrection *= -1.;
}

// Obtain iterator over the Correction per curve segment.
const Correction* getPathCorrections(MatrixWorkspace& data)
{
    static_assert(
        sizeof(Correction) == sizeof(Real) * GEODESIC_DOF,
        "Invalid size of geodesic correction.");
    if (data.pathCorrection.size() * sizeof(Real) !=
        data.nCurves * sizeof(Correction)) {
        throw std::runtime_error("Invalid size of path corrections vector.");
    }
    return reinterpret_cast<const Correction*>(&data.pathCorrection[0]);
}

} // namespace

//==============================================================================
//                      APPLYING GEODESIC CORRECTIONS
//==============================================================================

namespace
{

// Given a correction vector computed from minimizing the cost function,
// compute the maximum allowed stepsize along that correction vector.
// The allowed step is computed by approximating the surface locally as
// a circle in each direction, and limiting the radial displacement on that
// circle.
void calcMaxAllowedCorrectionStepSize(
    const CurveSegment& curve,
    const State& s,
    const Correction& c,
    Real maxAngularDisplacementInDegrees,
    Real& maxAllowedStepSize)
{
    auto UpdateMaxStepSize = [&](Real maxDisplacementEstimate, Real curvature) {
        const Real maxAngle = maxAngularDisplacementInDegrees / 180. * M_PI;
        const Real maxAllowedDisplacement = maxAngle / curvature;
        const Real allowedStepSize =
            std::abs(maxAllowedDisplacement / maxDisplacementEstimate);
        maxAllowedStepSize = std::min(maxAllowedStepSize, allowedStepSize);
    };

    const InstanceEntry& cache = curve.getInstanceEntry(s);

    // Clamp tangential displacement at the initial contact point.
    {
        const Real dxEst = c[0];
        const Real k     = cache.curvatures_P[0];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp binormal displacement at the initial contact point.
    {
        const Real dxEst = c[1];
        const Real k     = cache.curvatures_P[1];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp tangential displacement at the final contact point.
    {
        const Real dxEst = std::abs(c[0]) + std::abs(c[3]);
        const Real k     = cache.curvatures_Q[0];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp binormal displacement at the final contact point.
    {
        const Real a     = cache.jacobi_Q[0];
        const Real r     = cache.jacobi_Q[1];
        const Real dxEst = std::abs(c[1] * a) + std::abs(c[2] * r);
        const Real k     = cache.curvatures_Q[1];
        UpdateMaxStepSize(dxEst, k);
    }
}

// Apply the correction to the initial condition of the geodesic, and
// shoot a new geodesic, updating the cache variable.
void applyGeodesicCorrection(
    const CurveSegment& curve,
    const State& s,
    const Correction& c,
    const IntegratorTolerances& tols)
{
    // Get the previous geodesic.
    const InstanceEntry& cache = curve.getInstanceEntry(s);
    const FrenetFrame& X_SP    = cache.X_SP;

    const Real tau   = cache.torsion_P;
    const Real kappa = cache.curvatures_P[0];

    const UnitVec3& t = getTangent(X_SP);
    const UnitVec3& n = getNormal(X_SP);
    const UnitVec3& b = getBinormal(X_SP);

    // Get corrected initial conditions.
    const Vec3 dx               = t * c[0] + b * c[1];
    const Vec3 correctedPoint_S = X_SP.p() + dx;

    const Vec3 w                  = (kappa * c[0] - tau * c[1]) * b - c[2] * n;
    const Vec3 correctedTangent_S = t + cross(w, t);

    // Take the length correction, and add to the current length.
    const Real dl = c[3]; // Length increment is the last element.
    const Real correctedLength =
        std::max(cache.length + dl, 0.); // Clamp length to be nonnegative.

    // Shoot the new geodesic.
    curve.calcLocalGeodesic(
        s,
        correctedPoint_S,
        correctedTangent_S,
        correctedLength,
        cache.integratorInitialStepSize,
        tols);
}

} // namespace

//==============================================================================
//                      CABLE SPAN POSITION LEVEL CACHE
//==============================================================================

void CableSpan::Impl::calcPosInfo(const State& s, PosInfo& ppe) const
{
    // Path origin and termination point.
    const Vec3 x_O =
        getOriginBody().getBodyTransform(s).shiftFrameStationToBase(
            m_OriginPoint);
    const Vec3 x_I =
        getTerminationBody().getBodyTransform(s).shiftFrameStationToBase(
            m_TerminationPoint);

    ppe.originPoint_G      = x_O;
    ppe.terminationPoint_G = x_I;

    // Axes considered when computing the path error.
    const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

    ppe.loopIter = 0;
    while (true) {
        // Make sure all curve segments are realized to position stage.
        // This will transform all last computed geodesics to Ground frame, and
        // will update each curve's WrappingStatus.
        for (ObstacleIndex ix(0); ix < getNumCurveSegments(); ++ix) {
            getCurveSegment(ix).realizePosition(
                s,
                findPrevPoint(s, ix),
                findNextPoint(s, ix),
                getIntegratorTolerances());
        }

        // Now that the WrappingStatus of all curve segments is known: Count
        // the number of obstacles in contact with the path.
        const int nActive = countActive(s);

        // If the path contains no curved segments it is a straight line.
        if (nActive == 0) {
            // Update cache entry and stop solver.
            ppe.cableLength          = (x_I - x_O).norm();
            ppe.terminationTangent_G = ppe.originTangent_G =
                UnitVec3(x_I - x_O);
            break;
        }

        // If some obstacles are in contact with the cable the path error needs
        // to be checked. If the path error is small, i.e. there are no "kinks"
        // anywhere, the current path is OK. If the path error is too large,
        // corrections need to be computed for each curve segment in order to
        // drive the path error to zero.

        // Grab the shared data cache for helping with computing the path
        // corrections. This data is only used as an intermediate variable, and
        // will be discarded after each iteration. Note that the number active
        // segments determines the sizes of the matrices involved.
        MatrixWorkspace& data =
            getSubsystem().getImpl().updSolverData(s).updOrInsert(nActive);

        // Compute the straight-line segments of this cable span.
        calcLineSegments(*this, s, x_O, x_I, data.lineSegments);

        // Evaluate the current path error as the misalignment of the straight
        // line segments with the curve segment's tangent vectors at the
        // contact points.
        calcPathErrorVector<2>(s, data.lineSegments, axes, data.pathError);
        ppe.pathError = data.pathError.normInf();

        // Stop iterating if max path error is small, or max iterations has been
        // reached.
        if (ppe.pathError < m_PathAccuracy || ppe.loopIter >= m_SolverMaxIterations) {
            // Update cache entry and stop solver.
            ppe.cableLength =
                calcTotalCableLength(*this, s, data.lineSegments);
            ppe.originTangent_G      = data.lineSegments.front().direction;
            ppe.terminationTangent_G = data.lineSegments.back().direction;
            break;
        }

        // If the path error is too large corrections need to be applied to
        // each curve segment. Using the jacobian of the path error the
        // corrections can be computed that should drive the path error to
        // zero.

        // Evaluate the path error jacobian to the natural geodesic corrections
        // of each curve segment.
        calcPathErrorJacobian<2>(
            s,
            data.lineSegments,
            axes,
            data.pathErrorJacobian);

        // Compute the geodesic corrections for each curve segment: This gives
        // us a correction vector in a direction that lowers the path error.
        calcPathCorrections(data, ppe.pathError);

        // Compute the maximum allowed step size that we take along the
        // correction vector.
        Real stepSize            = 1.;
        const Correction* corrIt = getPathCorrections(data);
        forEachActiveCurveSegment(*this, s, [&](const CurveSegment& curve) {
            calcMaxAllowedCorrectionStepSize(
                curve,
                s,
                *corrIt,
                getMaxRadialStepInDegrees(),
                stepSize);
            ++corrIt;
        });
        data.pathCorrection *= stepSize;

        // Apply corrections to the curve segments. This will trigger shooting
        // new geodesics over the obstacles.
        corrIt = getPathCorrections(data);
        forEachActiveCurveSegment(*this, s, [&](const CurveSegment& curve) {
            applyGeodesicCorrection(
                curve,
                s,
                *corrIt,
                getIntegratorTolerances());
            ++corrIt;
        });

        // The applied corrections have changed the path: invalidate each
        // segment's cache.
        for (const CurveSegment& curve : m_CurveSegments) {
            // Also invalidate non-active segments: They might touchdown
            // again.
            curve.invalidatePosEntry(s);
        }

        ++ppe.loopIter;
    }
}

//==============================================================================
//                      CABLE SPAN VELOCITY LEVEL CACHE
//==============================================================================

void CableSpan::Impl::calcVelInfo(const State& s, VelInfo& velInfo) const
{
    const PosInfo& cablePE = getPosInfo(s);

    auto CalcPointVelocityInGround = [&](const MobilizedBody& mobod,
                                         const Vec3& point_G) -> Vec3 {
        // Not using MobilizedBody::findStationVelocityInGround because the
        // point is already in ground frame.

        // Get body kinematics in ground frame.
        const Vec3& x_BG = mobod.getBodyOriginLocation(s);
        const Vec3& w_BG = mobod.getBodyAngularVelocity(s);
        const Vec3& v_BG = mobod.getBodyOriginVelocity(s);

        // Compute surface point velocity in ground frame.
        return v_BG + w_BG % (point_G - x_BG);
    };

    Real& lengthDot = (velInfo.lengthDot = 0.);

    Vec3 v_GQ = getOriginBody().findStationVelocityInGround(s, m_OriginPoint);
    Vec3 x_GQ = m_OriginPoint;

    forEachActiveCurveSegment(*this, s, [&](const CurveSegment& curve) {
        const MobilizedBody& mobod            = curve.getMobilizedBody();
        const CurveSegment::PosEntry& curvePE = curve.getPosInfo(s);

        const UnitVec3& e_G = getTangent(curvePE.X_GP);

        const Vec3 v_GP = CalcPointVelocityInGround(mobod, curvePE.X_GP.p());

        lengthDot += dot(e_G, v_GP - v_GQ);

        x_GQ = curvePE.X_GQ.p();
        v_GQ = CalcPointVelocityInGround(mobod, x_GQ);
    });

    const Vec3 v_GP =
        getTerminationBody().findStationVelocityInGround(s, m_TerminationPoint);

    const UnitVec3 e_G(cablePE.terminationPoint_G - x_GQ);

    lengthDot += dot(e_G, v_GP - v_GQ);
}

//==============================================================================
//                               CURVE SEGMENT
//==============================================================================

CurveSegment::CurveSegment(
    CableSubsystem* subsystem,
    MobilizedBodyIndex body,
    const Transform& X_BS,
    ContactGeometry geometry,
    Vec3 initPointGuess) :
    m_Subsystem(subsystem),
    m_Body(body), m_X_BS(X_BS), m_Geometry(geometry),
    m_ContactPointHint_S(initPointGuess)
{
    SimTK_ASSERT(
        body.isValid(),
        "Failed to create new CurveSegment: Invalid MobilizedBodyIndex.");
}

const MobilizedBody& CurveSegment::getMobilizedBody() const
{
    return getSubsystem()
        .getImpl()
        .getMultibodySystem()
        .getMatterSubsystem()
        .getMobilizedBody(m_Body);
}

//==============================================================================
//                      SHOOTING GEODESICS HELPER
//==============================================================================

namespace
{

// Helper for shooting a new geodesic for this curve segment.
// The resulting geodesic is written to the provided cache variable.
void shootNewGeodesic(
    const CurveSegment& curve,
    Vec3 point_S,
    Vec3 tangent_S,
    Real length,
    Real initIntegratorStepSize,
    const IntegratorTolerances& tols,
    InstanceEntry& cache)
{
    cache.samples.clear();
    cache.integratorInitialStepSize = initIntegratorStepSize;

    const ContactGeometry& geometry = curve.getContactGeometry();

    // TODO implement shooting Analytic and Parametric geodesics.
    geometry.shootGeodesicInDirectionImplicitly(
        point_S,
        tangent_S,
        length,
        tols.intergatorAccuracy,
        tols.constraintProjectionTolerance,
        tols.constraintProjectionMaxIterations,
        cache.integratorInitialStepSize,
        cache.jacobi_Q,
        cache.jacobiDot_Q,
        // Provide function for logging the frenet frames during integration.
        [&](const Real& l, const Vec3& p, const Vec3& t) {
            // Compute frenet frame from position, normal and tangent.
            FrenetFrame X_Gk;
            X_Gk.setP(p);
            X_Gk.updR().setRotationFromTwoAxes(
                geometry.calcSurfaceUnitNormal(p),
                NormalAxis,
                t,
                TangentAxis);
            // Add frame to logged samples.
            cache.samples.emplace_back(l, X_Gk);
        });

    cache.X_SP = cache.samples.front().frame;
    cache.X_SQ = cache.samples.back().frame;

    cache.curvatures_P = {
        calcSurfaceCurvature(geometry, cache.X_SP, TangentAxis),
        calcSurfaceCurvature(geometry, cache.X_SP, BinormalAxis),
    };

    cache.curvatures_Q = {
        calcSurfaceCurvature(geometry, cache.X_SQ, TangentAxis),
        calcSurfaceCurvature(geometry, cache.X_SQ, BinormalAxis),
    };

    cache.torsion_P = geometry.calcSurfaceTorsionInDirection(
        cache.X_SP.p(),
        getTangent(cache.X_SP));
    cache.torsion_Q = geometry.calcSurfaceTorsionInDirection(
        cache.X_SQ.p(),
        getTangent(cache.X_SQ));

    cache.length = length;
}

} // namespace

//==============================================================================
//                  CURVE SEGMENT STATUS HANDLING
//==============================================================================
// Helpers for:
// - Shooting the initial geodesic.
// - Asserting path validity.
// - Detecting surface liftoff.
// - Detecting surface touchdown.

// Section containing helpers for curveSegment status handling.
namespace
{

// Helper function for computing the initial wrapping path as a zero length
// geodesic.
// TODO use the first four stages?
void calcInitCurveIfNeeded(
    const CurveSegment& curve,
    const State& s,
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S,
    const IntegratorTolerances& tols)
{
    if (curve.getInstanceEntry(s).status != WrappingStatus::InitialGuess) {
        return;
    }
    curve.calcLocalGeodesic(
        s,
        curve.getContactPointHint(),
        nextPoint_S - prevPoint_S,
        0.,
        NaN,
        tols);
}

// Assert that previous and next point lie above the surface. Points are in
// surface coordinates.
void assertSurfaceBounds(
    const CurveSegment& curve,
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S)
{
    // Try block was added here because some surfaces throw an exception when
    // calling calcSurfaceValue outside of a given range.
    auto IsPointAboveSurface = [&](const Vec3 point) -> bool {
        try {
            // TODO surface value is negative above surface.
            return curve.getContactGeometry().calcSurfaceValue(point) < 0.;
        } catch (...) {
            // TODO this is a bit too risky.
            return true;
        }
    };

    // Make sure that the previous point does not lie inside the surface.
    SimTK_ERRCHK_ALWAYS(
        IsPointAboveSurface(prevPoint_S),
        "CurveSegment::Impl::assertSurfaceBounds",
        "Preceding point lies inside the surface");
    SimTK_ERRCHK_ALWAYS(
        IsPointAboveSurface(nextPoint_S),
        "CurveSegment::Impl::assertSurfaceBounds",
        "Next point lies inside the surface");
}

void calcCurveTouchdownIfNeeded(
    const CurveSegment& curve,
    const State& s,
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S,
    const IntegratorTolerances& tols)
{
    // Given the straight line between the point before and after this segment,
    // touchdown is detected by computing the point on that line that is
    // nearest to the surface.

    // Only attempt touchdown when lifted.
    const InstanceEntry& cache = curve.getInstanceEntry(s);
    if (cache.status != WrappingStatus::LiftedFromSurface) {
        return;
    }

    // Use the cached tracking point as the initial guess.
    Vec3 pointOnLineNearSurface_S = cache.trackingPointOnLine_S;

    // Compute the point on the line nearest the surface.
    const bool touchdownDetected =
        curve.getContactGeometry().calcNearestPointOnLineImplicitly(
            prevPoint_S,
            nextPoint_S,
            tols.constraintProjectionMaxIterations,
            tols.constraintProjectionTolerance,
            pointOnLineNearSurface_S);

    // In case of touchdown, shoot a zero-length geodesic at the touchdown
    // point.
    if (touchdownDetected) {
        curve.calcLocalGeodesic(
            s,
            pointOnLineNearSurface_S,
            nextPoint_S - prevPoint_S,
            0.,
            cache.integratorInitialStepSize,
            tols);
        return;
    }

    // Not touchingdown indicates liftoff:
    curve.liftCurveFromSurface(s, pointOnLineNearSurface_S);
}

void calcCurveLiftoffIfNeeded(
    const CurveSegment& curve,
    const State& s,
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S)
{
    // Only attempt liftoff when currently wrapping the surface.
    const InstanceEntry& cache = curve.getInstanceEntry(s);
    if (cache.status != WrappingStatus::InContactWithSurface) {
        return;
    }

    // The curve length must have shrunk completely before lifting off.
    if (cache.length > 0.) {
        return;
    }

    // For a zero-length curve, trigger liftoff when the prev and next points
    // lie above the surface plane.
    if (dot(prevPoint_S - cache.X_SP.p(),
            cache.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0. ||
        dot(nextPoint_S - cache.X_SP.p(),
            cache.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0.) {
        // No liftoff.
        return;
    }

    // Liftoff detected, initialize the tracking point from the last contact
    // point.
    curve.liftCurveFromSurface(s, cache.X_SP.p());
}

} // namespace

void CurveSegment::calcLocalGeodesic(
    const State& s,
    Vec3 point_S,
    Vec3 tangent_S,
    Real length,
    Real stepSizeHint,
    const IntegratorTolerances& tols) const
{
    InstanceEntry& cache = updInstanceEntry(s);

    cache.status = WrappingStatus::InContactWithSurface;

    shootNewGeodesic(
        *this,
        point_S,
        tangent_S,
        length,
        stepSizeHint,
        tols,
        cache);

    getSubsystem().markDiscreteVarUpdateValueRealized(s, m_InstanceIx);

    invalidatePosEntry(s);
}

void CurveSegment::liftCurveFromSurface(const State& s, Vec3 trackingPoint_S)
    const
{
    InstanceEntry& cache = updInstanceEntry(s);

    cache.status                = WrappingStatus::LiftedFromSurface;
    cache.trackingPointOnLine_S = trackingPoint_S;

    getSubsystem().markDiscreteVarUpdateValueRealized(s, m_InstanceIx);

    invalidatePosEntry(s);
}

//==============================================================================
//                  CURVE SEGMENT: REALIZE POSITION LEVEL CACHE
//==============================================================================

void CurveSegment::realizePosition(
    const State& s,
    Vec3 prevPoint_G,
    Vec3 nextPoint_G,
    const IntegratorTolerances& tols) const
{
    if (getSubsystem().isCacheValueRealized(s, m_PosIx)) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "expected not realized when calling realizePosition");
    }

    if (getInstanceEntry(s).status == WrappingStatus::Disabled) {
        return;
    }

    // Compute tramsform from local surface frame to ground.
    const Transform X_GS = calcSurfaceFrameInGround(s);

    // Given the previous and next point determine if the cable wraps over
    // the obstacle surface at all.
    {
        // Transform the prev and next path points to the surface frame.
        const Vec3 prevPoint_S = X_GS.shiftBaseStationToFrame(prevPoint_G);
        const Vec3 nextPoint_S = X_GS.shiftBaseStationToFrame(nextPoint_G);

        // Make sure that the previous and next point do not lie inside the
        // surface body.
        assertSurfaceBounds(*this, prevPoint_S, nextPoint_S);

        calcInitCurveIfNeeded(*this, s, prevPoint_S, nextPoint_S, tols);

        // If previously not in contact with the surface, determine if we come
        // in contact now.
        calcCurveTouchdownIfNeeded(*this, s, prevPoint_S, nextPoint_S, tols);

        // When in contact with the surface, determine if we will lift off from
        // the surface now.
        calcCurveLiftoffIfNeeded(*this, s, prevPoint_S, nextPoint_S);
    }

    // At this point we have a valid geodesic in surface frame.
    const InstanceEntry& ie = getInstanceEntry(s);

    // Start updating the position level cache.
    PosEntry& ppe = updPosInfo(s);

    // Transform geodesic in local surface coordinates to ground.
    {
        // Store the local geodesic in ground frame.
        ppe.X_GS = X_GS;

        // Store the local geodesic in ground frame.
        ppe.X_GP = X_GS.compose(ie.X_SP);
        ppe.X_GQ = X_GS.compose(ie.X_SQ);
    }

    getSubsystem().markCacheValueRealized(s, m_PosIx);
}

//==============================================================================
//                      RESAMPLING GEODESIC PATH POINTS
//==============================================================================

// Helper functions for resampling the frenet frames obtained from shooting the
// geodesic.
namespace
{

// Compute the y at x using cubic Hermite interpolation between (y0, y0Dot) at
// x0, and (y1, y1Dot) at x1.
Vec3 calcHermiteInterpolation(
    Real x0,
    const Vec3& y0,
    const Vec3& y0Dot,
    Real x1,
    const Vec3& y1,
    const Vec3& y1Dot,
    Real x)
{
    // Interpolation is done by computing the coefficients of the polynomial:
    // f(s) = c0 + c1 * s + c2 * s^2 + c3 * s^3
    // followed by evaluating the polynomial for s = x - x0.

    const Real dx    = x1 - x0;
    const Vec3 dy    = y1 - y0;
    const Vec3 dyDot = y1Dot - y0Dot;

    // Compute the coefficients of the polynomial.
    const Vec3& c0 = y0;
    const Vec3& c1 = y0Dot;
    const Vec3 c3 =
        -2. * (dy - y0Dot * dx - 0.5 * dx * dyDot) / std::pow(dx, 3);
    const Vec3 c2 = (dyDot / dx - 3. * c3 * dx) / 2.;

    // Evaluate the polynomial.
    const Real s = x - x0;
    return c0 + s * (c1 + s * (c2 + s * c3));
}

// Compute the point in between two geodesic samples using Hermite
// interpolation.
Vec3 calcInterpolatedPoint(
    const LocalGeodesicSample& a,
    const LocalGeodesicSample& b,
    Real l)
{
    return calcHermiteInterpolation(
        a.length,
        a.frame.p(),
        getTangent(a.frame),
        b.length,
        b.frame.p(),
        getTangent(b.frame),
        l);
}

// Compute the tangent in between two geodesic samples using Hermite
// interpolation.
Vec3 calcInterpolatedTangent(
    const ContactGeometry& geometry,
    const LocalGeodesicSample& a,
    const LocalGeodesicSample& b,
    Real l)
{
    // Helper for calculating the derivative of the tangent.
    auto CalcTangentDot = [&](const FrenetFrame& X_S) -> Vec3 {
        const Real k = calcSurfaceCurvature(geometry, X_S, TangentAxis);
        return k * getNormal(X_S);
    };

    return calcHermiteInterpolation(
        a.length,
        getTangent(a.frame),
        CalcTangentDot(a.frame),
        b.length,
        getTangent(b.frame),
        CalcTangentDot(b.frame),
        l);
}

using InterpolatorFn = std::function<
    void(const LocalGeodesicSample& a, const LocalGeodesicSample& b, Real l)>;

// Resamples the geodesic samples at equal length intervals.
// @param geodesic contains the original geodesic samples.
// @param nSamples the number of samples to obtain.
// @param interpolator function for performing the interpolation between two
// selected samples at a given length.
// @return the number of samples used for the interpolation.
int resampleGeodesic(
    const std::vector<LocalGeodesicSample> geodesic,
    int nSamples, // TODO use length increment?
    InterpolatorFn& interpolator)
{
    // Some sanity checks.
    if (geodesic.empty()) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Resampling of geodesic failed: Provided geodesic is empty.");
    }
    if (geodesic.front().length != 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: First frame "
                                 "must be at length = zero");
    }
    if (geodesic.back().length < 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: Last frame "
                                 "must be at length > zero");
    }
    if (nSamples == 0) {
        throw std::runtime_error(
            "Invalid parameter: nSamples must be larger than zero.");
    }

    if (nSamples == 1 && geodesic.size() != 1) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: Requested "
                                 "number of samples must be unequal to 1");
    }

    // If there is but one sample in the geodesic, write that sample and exit.
    if (geodesic.size() == 1) {
        return 1;
    }

    // Seperate the interpolation points by equal length increments.
    const Real dl = geodesic.back().length / static_cast<Real>(nSamples - 1);

    // Compute the interpolated points from the geodesic.
    auto itGeodesic = geodesic.begin();
    // We can skip the first and last samples, because these are pushed
    // manually before and after this loop respectively (we start at i=1
    // and stop at i < nSamples-1).
    for (int i = 0; i < nSamples - 1; ++i) {

        // Length at the current interpolation point.
        const Real length = dl * static_cast<Real>(i);

        // Find the two samples (lhs, rhs) of the geodesic such that the
        // length of the interpolation point lies between them.
        // i.e. find: lhs.length <= length < rhs.length
        while (true) {
            // Sanity check: We should stay within range.
            if ((itGeodesic + 1) == geodesic.end()) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Attempted to read out of "
                    "array range");
            }

            // The candidate samples to use for interpolation.
            const LocalGeodesicSample& lhs = *itGeodesic;
            const LocalGeodesicSample& rhs = *(itGeodesic + 1);

            // Sanity check: Samples are assumed to be monotonically increasing
            // in length.
            if (lhs.length > rhs.length) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Samples are not "
                    "monotonically increasing in length.");
            }

            // Check that the interpolation point lies between these samples:
            // lhs.length <= length < rhs.length
            if (length >= rhs.length) {
                // Try the next two samples.
                ++itGeodesic;
                continue;
            }

            // Write interpolated point to the output buffer.
            interpolator(lhs, rhs, length);

            break;
        }
    }

    // Capture the last point of the geodesic.
    interpolator(geodesic.front(), geodesic.back(), geodesic.back().length);

    return nSamples;
}
} // namespace

int CurveSegment::calcPathPointsAndTangents(
    const State& state,
    const std::function<void(Real length, Vec3 point_G, UnitVec3 tangent_G)>& sink,
    int nSamples) const
{
    const InstanceEntry& geodesic_S = getInstanceEntry(state);
    if (!geodesic_S.isInContactWithSurface()) {
        return 0;
    }

    const CurveSegment::PosEntry& ppe = getPosInfo(state);

    // If the curve has zero length, return a single point, nothing to
    // interpolate.
    if (geodesic_S.length == 0.) {
        sink(0., ppe.X_GP.p(), getTangent(ppe.X_GP));
        return 1;
    }

    // Requesting two inerpolation points means we write the first and last
    // contact point. Nothing to interpolate.
    if (nSamples == 2) {
        sink(0., ppe.X_GP.p(), getTangent(ppe.X_GP));
        sink(geodesic_S.length, ppe.X_GQ.p(), getTangent(ppe.X_GQ));
        return 2;
    }

    // Resample the points from the integrator by interpolating at equal length
    // intervals.
    InterpolatorFn Interpolator = [&](const LocalGeodesicSample& a,
                                      const LocalGeodesicSample& b,
                                      Real l) {
        const Transform& X_GS = ppe.X_GS;
        const Vec3 interpolatedPoint_G =
            X_GS.shiftFrameStationToBase(calcInterpolatedPoint(a, b, l));
        const UnitVec3 interpolatedTangent_G =
            UnitVec3(X_GS.xformFrameVecToBase(
                calcInterpolatedTangent(getContactGeometry(), a, b, l)));

        sink(l, interpolatedPoint_G, interpolatedTangent_G);
    };
    return resampleGeodesic(geodesic_S.samples, nSamples, Interpolator);
}

int CurveSegment::calcPathPoints(
    const State& state,
    const std::function<void(Real length, Vec3 point_G)>& sink,
    int nSamples) const
{
    const InstanceEntry& geodesic_S = getInstanceEntry(state);
    if (!geodesic_S.isInContactWithSurface()) {
        return 0;
    }

    const CurveSegment::PosEntry& ppe = getPosInfo(state);

    // If the curve has zero length, return a single point, nothing to
    // interpolate.
    if (geodesic_S.length == 0.) {
        sink(0., ppe.X_GP.p());
        return 1;
    }

    // Requesting two inerpolation points means we write the first and last
    // contact point. Nothing to interpolate.
    if (nSamples == 2) {
        sink(0., ppe.X_GP.p());
        sink(geodesic_S.length, ppe.X_GQ.p());
        return 2;
    }

    // Resample the points from the integrator by interpolating at equal
    // intervals.

    // Define the required function that interpolates between two samples, and
    // logs the interpolated sample.
    InterpolatorFn Interpolator = [&](const LocalGeodesicSample& a,
                                      const LocalGeodesicSample& b,
                                      Real l) {
        const Vec3 interpolatedPoint_S = calcInterpolatedPoint(a, b, l);
        const Vec3 interpolatedPoint_G =
            ppe.X_GS.shiftFrameStationToBase(interpolatedPoint_S);
        sink(l, interpolatedPoint_G);
    };
    return resampleGeodesic(geodesic_S.samples, nSamples, Interpolator);
}

//==============================================================================
//                      SUBSYSTEM TESTING HELPER
//==============================================================================

bool CableSubsystemTestHelper::applyPerturbationTest(
    const MultibodySystem& system,
    const CableSubsystem& subsystem,
    const State& s,
    Real perturbation,
    Real bound,
    std::ostream& os)
{
    // Output of this perturbation test.
    bool success = true;

    for (const CableSpan& cableIt : subsystem.getImpl().cables) {
        const CableSpan::Impl& cable = cableIt.getImpl();

        // Axes considered when computing the path error.
        const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

        // We do not want to mess with the actual state, so we make a copy.
        const State sCopy = s;
        system.realize(sCopy, Stage::Position);

        // Count the number of active curve segments.
        const int nActive = cable.countActive(sCopy);

        if (nActive == 0) {
            os << "Cable has no active segments: Skipping perturbation test\n";
            continue;
        }

        MatrixWorkspace& data =
            subsystem.getImpl().updSolverData(sCopy).updOrInsert(nActive);

        constexpr int DOF = 4;
        for (int i = 0; i < (DOF * 2 * nActive); ++i) {
            cable.invalidatePositionLevelCache(sCopy);

            // Trigger realizing position level cache, resetting the
            // configuration.
            const CableSpan::Impl::PosInfo& ppe = cable.getPosInfo(sCopy);

            // Define the perturbation we will use for testing the jacobian.
            perturbation *= -1.;
            data.pathCorrection.setTo(0.);
            data.pathCorrection.set(i / 2, perturbation);

            calcLineSegments(
                cable,
                sCopy,
                ppe.originPoint_G,
                ppe.terminationPoint_G,
                data.lineSegments);

            cable.calcPathErrorVector<2>(
                sCopy,
                data.lineSegments,
                axes,
                data.pathError);

            cable.calcPathErrorJacobian<2>(
                sCopy,
                data.lineSegments,
                axes,
                data.pathErrorJacobian);

            Vector predictedPathError =
                data.pathErrorJacobian * data.pathCorrection + data.pathError;

            std::vector<WrappingStatus> prevWrappingStatus{};
            for (const CurveSegment& curve : cable.m_CurveSegments) {
                prevWrappingStatus.push_back(curve.getInstanceEntry(sCopy).status);
            }

            const Correction* corrIt = getPathCorrections(data);
            forEachActiveCurveSegment(
                cable,
                sCopy,
                [&](const CurveSegment& curve) {
                    applyGeodesicCorrection(
                        curve,
                        sCopy,
                        *corrIt,
                        cable.getIntegratorTolerances());
                    ++corrIt;
                });

            for (const CurveSegment& curve : cable.m_CurveSegments) {
                curve.invalidatePosEntry(sCopy);
            }

            for (ObstacleIndex ix(0); ix < cable.getNumCurveSegments(); ++ix) {
                cable.getCurveSegment(ix).realizePosition(
                    sCopy,
                    cable.findPrevPoint(sCopy, ix),
                    cable.findNextPoint(sCopy, ix),
                    cable.getIntegratorTolerances());
            }

            calcLineSegments(
                cable,
                sCopy,
                ppe.originPoint_G,
                ppe.terminationPoint_G,
                data.lineSegments);
            cable.calcPathErrorVector<2>(
                sCopy,
                data.lineSegments,
                axes,
                data.pathError);

            bool WrappingStatusChanged = false;
            {
                int ix = -1;
                for (const CurveSegment& curve : cable.m_CurveSegments) {
                    WrappingStatusChanged |= prevWrappingStatus.at(++ix) !=
                                             curve.getInstanceEntry(sCopy).status;
                }
            }

            if (WrappingStatusChanged) {
                os << "Wrapping status changed: Stopping test\n";
                break;
            }

            const Real predictionError =
                (predictedPathError - data.pathError).norm();

            bool passedTest = predictionError / perturbation <= bound;
            if (!passedTest) {
                os << "FAILED perturbation test for correction = "
                   << data.pathCorrection << "\n";
            } else {
                os << "PASSED perturbation test for correction = "
                   << data.pathCorrection << "\n";
            }
            os << "    Got      : " << data.pathError << "\n";
            os << "    Predicted: " << predictedPathError << "\n";
            os << "    Err      : "
               << (predictedPathError - data.pathError) / perturbation << "\n";
            os << "    Err norm : "
               << (predictedPathError - data.pathError).normInf() / perturbation
               << "\n";
            success &= passedTest;
        }
    }
    return success;
}
