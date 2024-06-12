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

#include "simbody/internal/CableSpan.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

using namespace SimTK;

using WrappingStatus = CableSpan::ObstacleWrappingStatus;
using ObstacleIndex  = CableSpan::ObstacleIndex;

//==============================================================================
//                          Helper Structures
//==============================================================================
namespace
{

//==============================================================================
//                          Struct Line Segment
//==============================================================================
// Representation of a cable segment that does not lie on a surface: A straight
// line.
struct LineSegment final
{
    LineSegment() = default;

    // Construct a new LineSegment connecting pointA and pointB.
    LineSegment(const Vec3& pointA, const Vec3& pointB) :
        length((pointB - pointA).norm()), direction((pointB - pointA) / length)
    {}

    Real length = NaN;
    UnitVec3 direction{NaN, NaN, NaN};
};

//==============================================================================
//                         Struct MatrixWorkspace
//==============================================================================
/** This is a helper struct that is used by a CableSpan to compute the
Stage::Position level data, i.e. the spanned path.
Computing the spanned path involves a Newton type iteration, for which several
matrices need to be computed. This struct provides the required matrices.
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
//                            Frenet Frame
//==============================================================================
// The frenet frame plays an important role in defining the state of a geodesic.
// In the end it is simply a Transform, which is why we define it as an alias.
//
// The position of the frenet frame is defined to lie on the geodesic, but for
// the orientation of the frame different conventions are used.
// Here we define it as:
// - X axis: tangent to geodesic
// - Y axis: normal to surface
// - Z axis: binormal to geodesic (tangent cross normal)
using FrenetFrame                        = Transform;
static const CoordinateAxis TangentAxis  = CoordinateAxis::XCoordinateAxis();
static const CoordinateAxis NormalAxis   = CoordinateAxis::YCoordinateAxis();
static const CoordinateAxis BinormalAxis = CoordinateAxis::ZCoordinateAxis();

// Helper for grabbing the tangent axis.
const UnitVec3& getTangent(const FrenetFrame& X)
{
    return X.R().getAxisUnitVec(TangentAxis);
}

// Helper for grabbing the normal axis.
const UnitVec3& getNormal(const FrenetFrame& X)
{
    return X.R().getAxisUnitVec(NormalAxis);
}

// Helper for grabbing the binormal axis.
const UnitVec3& getBinormal(const FrenetFrame& X)
{
    return X.R().getAxisUnitVec(BinormalAxis);
}

// Helper for computing the frenet frame at a knot point of a geodesic.
FrenetFrame calcFrenetFrameFromGeodesicState(
    const ContactGeometry& geometry,
    const ContactGeometry::ImplicitGeodesicState& q)
{
    FrenetFrame X;
    X.updP() = q.point;
    X.updR().setRotationFromTwoAxes(
        geometry.calcSurfaceUnitNormal(q.point),
        NormalAxis,
        q.tangent,
        TangentAxis);
    return X;
}

} // namespace

//==============================================================================
//                      Cached Data Structures
//==============================================================================
// The following section contains data structures that are cached during a
// simulation.
namespace
{

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
    // Constructs requested MatrixWorkspace of requested dimension if not done
    // previously.
    MatrixWorkspace& updOrInsert(int nActive)
    {
        SimTK_ASSERT1(
            nActive > 0,
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

namespace CableSpanData
{
//=============================================================================
//                          CableSpanData::Pos
//=============================================================================
/** CurveSegmentData::Pos is a data structure used by class CableSpan::Impl to
store quantities that can be computed at Stage::Position.

Member variables with an "_G" suffix are expressed in the Ground frame. */
struct Pos final
{
    // Position of the cable origin point in the ground frame.
    Vec3 originPoint_G{NaN, NaN, NaN};
    // Position of the cable termination point in the ground frame.
    Vec3 terminationPoint_G{NaN, NaN, NaN};
    // Tangent vector to the cable at the origin point, in the ground frame.
    UnitVec3 originTangent_G{NaN, NaN, NaN};
    // Tangent vector to the cable at the termination point, in the ground
    // frame.
    UnitVec3 terminationTangent_G{NaN, NaN, NaN};
    // The total cable length.
    Real cableLength = NaN;
    // The maximum path error.
    Real pathError = NaN;
    // Number of iterations the solver used to compute the cable's path.
    int loopIter = -1;
};

//=============================================================================
//                          CableSpanData::Vel
//=============================================================================
/** CurveSegmentData::Vel is a data structure used by class CableSpan::Impl to
store quantities that can be computed at Stage::Velocity.
TODO currently there is not much in it. */
struct Vel final
{
    // Time derivative of the total cable length.
    Real lengthDot = NaN;
};

} // namespace CableSpanData

namespace CurveSegmentData
{
/** CurveSegmentData::Instance is a data structure used by class CurveSegment to
store quantities that can be computed at Stage::Instance.

At the heart of computing the cable path lies a Newton iteration, which
greatly benefits from a good solution estimate. Using a discrete cache
entry allows for starting the solver from the previous solution.
Furthermore, by storing the previous curve segment, this allows winding over
obstacles multiple times. Therefore this struct contains information related to
the local (body fixed) geodesic.

A discrete auto update cache variable is used because the cable path over an
obstacle can be seen as a nonholonomic constraint. If we wish to redo
the cable path computation for some reason, we must be able to reset to
the previous path.

This data becomes available after Stage::Instance.

Following Scholz2015 the following subscripts are used:

Member variables with an "_P" suffix indicate the start of the geodesic.
Member variables with an "_Q" suffix indicate the end of the geodesic.
Member variables with an "_S" suffix are expressed in the Surface frame.

For example: point_SP would indicate the initial contact point represented and
measured from the surface's origin frame. */
struct Instance final
{
    // Frenet frame at the initial contact point w.r.t. the surface frame.
    FrenetFrame X_SP{};
    // Frenet frame at the final contact point w.r.t. the surface frame.
    FrenetFrame X_SQ{};

    // Length of this curve segment.
    Real length = NaN;

    // Curvatures at the contact points, with the tangential and binormal
    // direction as the first and last element respectively.
    Vec2 curvatures_P{NaN}; // initial contact point.
    Vec2 curvatures_Q{NaN}; // final contact point.

    // Geodesic torsion at the initial and final contact points.
    Real torsion_P = NaN;
    Real torsion_Q = NaN;

    // Jacobi scalars at the final contact point. Contains the
    // translational and rotational direction as the first and last
    // element (see Scholz2015).
    Vec2 jacobi_Q{NaN};
    Vec2 jacobiDot_Q{NaN};

    // Samples recorded from shooting the geodesic over the surface.
    // For an analytic contact geometry this container will be empty.
    // Otherwise, the first and last sample will contain X_SP and X_SQ
    // defined above.
    std::vector<ContactGeometry::ImplicitGeodesicState>
        geodesicIntegratorStates;

    // The initial integrator stepsize to try next time when shooting a
    // geodesic. This step size estimate will improve each time after
    // shooting a new geodesic.
    Real integratorInitialStepSize = NaN;

    // Given the line spanned by the point before and after this curve
    // segment, the tracking point lies on that line and is nearest to the
    // surface. This point is used to find the touchdown location when the
    // curve is not in contact with the obstacles's surface.
    Vec3 trackingPointOnLine_S{NaN, NaN, NaN};

    // This is a flag to indicate whether the cable comes into contact with
    // the obstacle at all. If the cable is not in contact with the
    // surface, only trackingPointOnLine_S will contain valid data. If the
    // cable is in contact with the surface, trackingPointOnLine_S will not
    // contain valid data.
    WrappingStatus status = WrappingStatus::InitialGuess;
};

// Stage::Position level cache entry for holding the curve information in
// ground frame.
struct Pos final
{
    const Vec3& getFirstContactPoint() const
    {
        return X_GP.p();
    }

    const Vec3& getFinalContactPoint() const
    {
        return X_GQ.p();
    }

    // Position and orientation of contact geometry w.r.t. ground.
    Transform X_GS{};

    // Frenet frame at the initial contact point w.r.t. ground.
    FrenetFrame X_GP{};
    // Frenet frame at the final contact point w.r.t. ground.
    FrenetFrame X_GQ{};
};

} // namespace CurveSegmentData

} // namespace

//==============================================================================
//                          Cable Span Parameters
//==============================================================================
namespace
{
// Helper struct for collecting all parameters for setting up the
// integrator that computes the geodesic over the obstacle's surface.
// TODO collect all parameters.
struct IntegratorTolerances
{
    Real intergatorAccuracy = 1e-6; // TODO set to reasonable
    // TODO combine constraintProjectionTolerance with intergatorAccuracy?
    Real constraintProjectionTolerance = 1e-6; // TODO set to reasonable
    // TODO this is not currently connected to anything...
    int constraintProjectionMaxIterations = 50; // TODO set to reasonable
};
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
        X_S.R().getAxisUnitVec(direction));
}

// Helper function for flipping the sign such that positive value is outside of
// the surface.
// TODO is not used?
Real calcSurfaceValue(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceValue(point_S);
}
// Flipping sign...
// TODO is not used?
Vec3 calcSurfaceGradient(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceGradient(point_S);
}
// Flipping sign...
// TODO is not used?
Mat33 calcSurfaceHessian(const ContactGeometry& geometry, const Vec3& point_S)
{
    return -geometry.calcSurfaceHessian(point_S);
}

} // namespace

//==============================================================================
//                      Class CableSubsystem::Impl
//==============================================================================
// TODO DESCRIPTION
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

    CableSpanIndex adoptCable(CableSpan& cable)
    {
        invalidateSubsystemTopologyCache();

        CableSpanIndex cableIx(cables.size());
        cables.emplace_back(cable);
        return cableIx;
    }

    int getNumCables() const
    {
        return cables.size();
    }

    const CableSpan& getCable(CableSpanIndex index) const;

    CableSpan& updCable(CableSpanIndex index);

    // Return the MultibodySystem which owns this WrappingPathSubsystem.
    const MultibodySystem& getMultibodySystem() const
    {
        return MultibodySystem::downcast(getSystem());
    }

    // TODO grab correct dimension here.
    PathSolverScratchData& updSolverData(const State& state) const
    {
        return Value<PathSolverScratchData>::updDowncast(
            updCacheEntry(state, m_CacheIx));
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    Impl* cloneImpl() const override
    {
        return new Impl(*this);
    }

    int realizeSubsystemTopologyImpl(State& state) const override;

    int calcDecorativeGeometryAndAppendImpl(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const override;

    // TOPOLOGY STATE
    Array_<CableSpan, CableSpanIndex> cables;

    CacheEntryIndex m_CacheIx;

    friend CableSubsystemTestHelper;
};

//==============================================================================
//                          Struct Curve Segment
//==============================================================================
class CurveSegment final
{
public:
    //--------------------------------------------------------------------------
    // Constructors.
    //--------------------------------------------------------------------------

    CurveSegment() = default;

    CurveSegment(
        CableSubsystem* subsystem,
        MobilizedBodyIndex body,
        const Transform& X_BS,
        ContactGeometry geometry,
        Vec3 initPointGuess) :
        m_Subsystem(subsystem),
        m_Body(body), m_X_BS(X_BS), m_Geometry(geometry),
        m_ContactPointHint_S(initPointGuess)
    {}

    //--------------------------------------------------------------------------
    // Realizing and cache access.
    //--------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& s)
    {
        // Allocate an auto-update discrete variable for the last computed
        // geodesic.
        m_InstanceIx = updSubsystem().allocateAutoUpdateDiscreteVariable(
            s,
            Stage::Report,
            new Value<CurveSegmentData::Instance>(),
            Stage::Position);

        m_PosIx = updSubsystem().allocateCacheEntry(
            s,
            Stage::Position,
            Stage::Infinity,
            new Value<CurveSegmentData::Pos>());

        m_Decoration = getContactGeometry()
                           .createDecorativeGeometry()
                           .setColor(Orange)
                           .setOpacity(.75)
                           .setResolution(3);
    }

    void realizePosition(
        const State& s,
        Vec3 prevPoint_G,
        Vec3 nextPoint_G,
        const IntegratorTolerances& tols) const;

    void invalidatePosEntry(const State& state) const
    {
        getSubsystem().markCacheValueNotRealized(state, m_PosIx);
    }

    //------------------------------------------------------------------------------

    // Set the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    void setContactPointHint(Vec3 contactPointHint_S)
    {
        m_ContactPointHint_S = contactPointHint_S;
    }

    // Get the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    Vec3 getContactPointHint() const
    {
        return m_ContactPointHint_S;
    }

    const ContactGeometry& getContactGeometry() const
    {
        return m_Geometry;
    }
    void setContactGeometry(ContactGeometry geometry)
    {
        m_Geometry = geometry;
    }

    const DecorativeGeometry& getDecoration() const
    {
        return m_Decoration;
    }

    const MobilizedBody& getMobilizedBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_Body);
    }

    const MobilizedBodyIndex& getMobilizedBodyIndex() const
    {
        return m_Body;
    }

    void setMobilizedBodyIndex(MobilizedBodyIndex body)
    {
        m_Body = body;
    }

    const Transform& getXformSurfaceToBody() const
    {
        return m_X_BS;
    }

    void setXformSurfaceToBody(const Transform& X_BS)
    {
        m_X_BS = X_BS;
    }

    const CableSubsystem& getSubsystem() const
    {
        if (!m_Subsystem) {
            SimTK_ERRCHK_ALWAYS(
                m_Subsystem,
                "CurveSegment::getSubsystem()",
                "CableSpan not yet adopted by any CableSubsystem");
        }
        return *m_Subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        if (!m_Subsystem) {
            SimTK_ERRCHK_ALWAYS(
                m_Subsystem,
                "CurveSegment::updSubsystem()",
                "CableSpan not yet adopted by any CableSubsystem");
        }
        return *m_Subsystem;
    }

    void setSubsystem(CableSubsystem& subsystem)
    {
        m_Subsystem = &subsystem;
    }

    //------------------------------------------------------------------------------

    const CurveSegmentData::Instance& getDataInst(const State& s) const
    {
        const CableSubsystem& subsystem = getSubsystem();
        if (!subsystem.isDiscreteVarUpdateValueRealized(s, m_InstanceIx)) {
            updDataInst(s) = getPrevInstanceEntry(s);
            subsystem.markDiscreteVarUpdateValueRealized(s, m_InstanceIx);
        }
        return Value<CurveSegmentData::Instance>::downcast(
            subsystem.getDiscreteVarUpdateValue(s, m_InstanceIx));
    }

    const CurveSegmentData::Pos& getDataPos(const State& s) const
    {
        return Value<CurveSegmentData::Pos>::downcast(
            getSubsystem().getCacheEntry(s, m_PosIx));
    }

    //------------------------------------------------------------------------------

    Transform calcSurfaceFrameInGround(const State& s) const
    {
        return getMobilizedBody().getBodyTransform(s).compose(m_X_BS);
    }

    bool isInContactWithSurface(const State& s) const
    {
        const CurveSegmentData::Instance& dataInst = getDataInst(s);
        return dataInst.status == WrappingStatus::InContactWithSurface ||
               dataInst.status == WrappingStatus::InitialGuess;
    }

    // See CurveSegment::calcPathPoints for description.
    int calcPathPoints(
        const State& state,
        const std::function<void(Real length, Vec3 point_G)>& sink,
        int nSamples = 0) const;

    int calcPathPointsAndTangents(
        const State& state,
        const std::function<
            void(Real length, Vec3 point_G, UnitVec3 tangent_G)>& sink,
        int nSamples) const;

    // Compute a new geodesic from provided initial conditions.
    // This method will update the Instance level cache, and invalidates the
    // Stage::Position level cache.
    void calcLocalGeodesic(
        const State& s,
        Vec3 point_S,
        Vec3 tangent_S,
        Real length,
        Real stepSizeHint,
        const IntegratorTolerances& tols) const;

    // Lift curve from surface, and start tracking the given point.
    // This method will update the Instance level cache, and invalidates the
    // Stage::Position level cache.
    void liftCurveFromSurface(const State& s, Vec3 trackingPoint_S) const;

    void storeCurrentPath(State& state) const
    {
        updPrevInstanceEntry(state) = getDataInst(state);
    }

    //------------------------------------------------------------------------------

private:
    CurveSegmentData::Instance& updDataInst(const State& state) const
    {
        return Value<CurveSegmentData::Instance>::updDowncast(
            getSubsystem().updDiscreteVarUpdateValue(state, m_InstanceIx));
    }

    const CurveSegmentData::Instance& getPrevInstanceEntry(const State& state) const
    {
        return Value<CurveSegmentData::Instance>::downcast(
            getSubsystem().getDiscreteVariable(state, m_InstanceIx));
    }

    CurveSegmentData::Instance& updPrevInstanceEntry(State& state) const
    {
        return Value<CurveSegmentData::Instance>::updDowncast(
            getSubsystem().updDiscreteVariable(state, m_InstanceIx));
    }

    CurveSegmentData::Pos& updPosInfo(const State& state) const
    {
        return Value<CurveSegmentData::Pos>::updDowncast(
            getSubsystem().updCacheEntry(state, m_PosIx));
    }

    //------------------------------------------------------------------------------

    // Subsystem info.
    CableSubsystem* m_Subsystem; // TODO manage reference count?

    // MobilizedBody that surface is attached to.
    MobilizedBodyIndex m_Body;
    // Surface to body transform.
    Transform m_X_BS;

    // Obstacle surface.
    ContactGeometry m_Geometry;
    DecorativeGeometry m_Decoration;

    // Topology cache.
    CacheEntryIndex m_PosIx;
    DiscreteVariableIndex m_InstanceIx;

    // Initial contact point hint used to setup the initial path.
    Vec3 m_ContactPointHint_S{NaN, NaN, NaN};

    // Helper class for unit tests.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                         Class Cablespan::Impl
//==============================================================================
class CableSpan::Impl
{
public:
    //--------------------------------------------------------------------------
    // Constructors.
    //--------------------------------------------------------------------------

    Impl() = default;

    // Construct a CableSpan::Impl with valid data fields, but without being
    // registered at any CableSubsystem.
    Impl(
        MobilizedBodyIndex originBody,
        Vec3 originPoint,
        MobilizedBodyIndex terminationBody,
        Vec3 terminationPoint) :
        m_OriginBody(originBody),
        m_OriginPoint(originPoint), m_TerminationBody(terminationBody),
        m_TerminationPoint(terminationPoint)
    {}

    //--------------------------------------------------------------------------
    // Realizing and accessing cache.
    //--------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& state)
    {
        indexDataPos = updSubsystem().allocateCacheEntry(
            state,
            Stage::Position,
            Stage::Infinity,
            new Value<CableSpanData::Pos>());

        indexDataVel = updSubsystem().allocateCacheEntry(
            state,
            Stage::Velocity,
            Stage::Infinity,
            new Value<CableSpanData::Vel>());

        for (CurveSegment& segment : m_CurveSegments) {
            segment.realizeTopology(state);
        }
    }

    void invalidateTopology()
    {
        if (m_Subsystem) {
            getSubsystem().invalidateSubsystemTopologyCache();
        }
    }

    void realizePosition(const State& state) const
    {
        if (getSubsystem().isCacheValueRealized(state, indexDataPos)) {
            return;
        }
        calcPosInfo(state, updPosInfo(state));
        getSubsystem().markCacheValueRealized(state, indexDataPos);
    }

    void realizeVelocity(const State& state) const
    {
        realizePosition(state);
        if (getSubsystem().isCacheValueRealized(state, indexDataVel)) {
            return;
        }
        calcVelInfo(state, updVelInfo(state));
        getSubsystem().markCacheValueRealized(state, indexDataVel);
    }

    const CableSpanData::Pos& getPosInfo(const State& state) const
    {
        realizePosition(state);
        return Value<CableSpanData::Pos>::downcast(
            getSubsystem().getCacheEntry(state, indexDataPos));
    }

    const CableSpanData::Vel& getVelInfo(const State& state) const
    {
        realizeVelocity(state);
        return Value<CableSpanData::Vel>::downcast(
            getSubsystem().getCacheEntry(state, indexDataVel));
    }

    //--------------------------------------------------------------------------
    // Obstacle access.
    //--------------------------------------------------------------------------

    ObstacleIndex addSurfaceObstacle(
        MobilizedBodyIndex mobod,
        const Transform& X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint)
    {
        invalidateTopology();

        ObstacleIndex obstacleIx(m_CurveSegments.size());

        m_CurveSegments.push_back(
            CurveSegment(m_Subsystem, mobod, X_BS, geometry, contactPointHint));

        return obstacleIx;
    }

    int getNumCurveSegments() const
    {
        return m_CurveSegments.size();
    }

    const CurveSegment& getCurveSegment(ObstacleIndex ix) const
    {
        return m_CurveSegments[ix];
    }

    CurveSegment& updCurveSegment(ObstacleIndex ix)
    {
        return m_CurveSegments[ix];
    }

    //--------------------------------------------------------------------------
    // Cable end points interface.
    //--------------------------------------------------------------------------

    // Get body to which this cable's origin point is rigidly attached to.
    const Mobod& getOriginBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_OriginBody);
    }

    // Get body to which this cable's termination point is rigidly attached to.
    const Mobod& getTerminationBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_TerminationBody);
    }

    // Get the origin point in body fixed coordinates.
    const Vec3 getOriginPoint_B() const
    {
        return m_OriginPoint;
    }

    // Get the termination point in body fixed coordinates.
    const Vec3 getTerminationPoint_B() const
    {
        return m_TerminationPoint;
    }

    //--------------------------------------------------------------------------
    // Cable parameters interface.
    //--------------------------------------------------------------------------

    // TODO collect params
    Real getSurfaceConstraintTolerance() const
    {
        return m_IntegratorTolerances.constraintProjectionTolerance;
    }
    void setSurfaceConstraintTolerance(Real tolerance)
    {
        m_IntegratorTolerances.constraintProjectionTolerance = tolerance;
    }

    int getSurfaceProjectionMaxIter() const
    {
        return m_IntegratorTolerances.constraintProjectionMaxIterations;
    }
    void setSurfaceProjectionMaxIter(int maxIterations)
    {
        m_IntegratorTolerances.constraintProjectionMaxIterations =
            maxIterations;
    }

    Real getIntegratorAccuracy() const
    {
        return m_IntegratorTolerances.intergatorAccuracy;
    }
    void setIntegratorAccuracy(Real accuracy)
    {
        m_IntegratorTolerances.intergatorAccuracy = accuracy;
    }

    int getSolverMaxIterations() const
    {
        return m_SolverMaxIterations;
    }
    void setSolverMaxIterations(int maxIterations)
    {
        m_SolverMaxIterations = maxIterations;
    }

    Real getPathAccuracy() const
    {
        return m_PathAccuracy;
    }
    void setPathAccuracy(Real accuracy)
    {
        m_PathAccuracy = accuracy;
    }

    Real getMaxRadialStepInDegrees() const
    {
        return m_MaxCorrectionStepDeg;
    }
    void setMaxRadialStepInDegrees(Real maxStepInDeg)
    {
        m_MaxCorrectionStepDeg = maxStepInDeg;
    }

    const IntegratorTolerances& getIntegratorTolerances() const
    {
        return m_IntegratorTolerances;
    }

    //--------------------------------------------------------------------------
    // CableSpan interface.
    //--------------------------------------------------------------------------

    // Find the number of CurveSegments that are in contact with an obstacle's
    // surface.
    int countActive(const State& s) const
    {
        int count = 0;
        for (const CurveSegment& segment : m_CurveSegments) {
            if (segment.isInContactWithSurface(s)) {
                ++count;
            }
        }
        return count;
    }

    // See CableSpan::calcPathPoints.
    int calcPathPoints(
        const State& state,
        Real lengthIncrement,
        const std::function<void(Real length, Vec3 point_G)>& sink) const
    {
        // Count number of points written.
        int count = 0;

        // Write the initial point.
        const CableSpanData::Pos& pos = getPosInfo(state);
        sink(0., pos.originPoint_G);
        ++count;

        // Write points along each of the curves.
        for (const CurveSegment& curve : m_CurveSegments) {

            Real curveLength = curve.getDataInst(state).length;
            int nPoints =
                static_cast<int>(std::abs(curveLength / lengthIncrement) + 1.);
            nPoints = nPoints < 2 ? 2 : nPoints;

            count += curve.calcPathPoints(state, sink, nPoints);
        }

        // Write the termination point.
        sink(pos.cableLength, pos.terminationPoint_G);
        ++count;

        // Return number of points written.
        return count;
    }

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    Real calcCablePower(const State& state, Real tension) const;

    int calcDecorativeGeometryAndAppend(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const;

private:
    // Get the subsystem this CableSpan is part of.
    const CableSubsystem& getSubsystem() const
    {
        SimTK_ERRCHK_ALWAYS(
            m_Subsystem,
            "CableSpan::Impl::getSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");
        return *m_Subsystem;
    }

    void setSubsystem(CableSubsystem& subsystem, CableSpanIndex index)
    {
        m_Subsystem = &subsystem;
        m_Index     = index;
        for (CurveSegment& curve : m_CurveSegments) {
            curve.setSubsystem(subsystem);
        }
    }

    CableSubsystem& updSubsystem()
    {
        SimTK_ERRCHK_ALWAYS(
            m_Subsystem,
            "CableSpan::Impl::updSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");
        return *m_Subsystem;
    }

    CableSpanIndex getIndex() const
    {
        SimTK_ERRCHK1_ALWAYS(
            m_Index.isValid(),
            "CableSpan::Impl::getIndex()",
            "Index %d is not valid.",
            m_Index);
        return m_Index;
    }

    CableSpanData::Pos& updPosInfo(const State& state) const
    {
        return Value<CableSpanData::Pos>::updDowncast(
            getSubsystem().updCacheEntry(state, indexDataPos));
    }

    CableSpanData::Vel& updVelInfo(const State& state) const
    {
        return Value<CableSpanData::Vel>::updDowncast(
            getSubsystem().updCacheEntry(state, indexDataVel));
    }

    void calcPosInfo(const State& s, CableSpanData::Pos& dataPos) const;
    void calcVelInfo(const State& s, CableSpanData::Vel& dataVel) const;

    //------------------------------------------------------------------------------
    //                      HELPER FUNCTIONS
    //------------------------------------------------------------------------------

    // Find the last contact point before the given curve segment, skipping
    // over any that are not in contact with their respective obstacle's
    // surface.
    Vec3 findPrevPoint(const State& state, ObstacleIndex ix) const
    {
        // Check if there is a curve segment preceding given obstacle.
        const CurveSegment* prevCurve = findPrevActiveCurveSegment(state, ix);
        if (prevCurve) {
            // If so, the previous point is the final contact point of the
            // previous curve.
            const Vec3& finalContactPoint_S =
                prevCurve->getDataInst(state).X_SQ.p();
            // Transform the contact point to ground.
            const Transform& X_GS = prevCurve->calcSurfaceFrameInGround(state);
            return X_GS.shiftFrameStationToBase(finalContactPoint_S);
        }
        // There are no curve segments before given obstacle: the previous point
        // is the path's origin point.
        return getOriginBody().getBodyTransform(state).shiftFrameStationToBase(
            m_OriginPoint);
    }

    // TODO DESCRIPTION
    Vec3 findNextPoint(const State& state, ObstacleIndex ix) const
    {
        // Check if there is a curve segment after given obstacle.
        const CurveSegment* nextCurve = findNextActiveCurveSegment(state, ix);
        if (nextCurve) {
            // If so, the next point is the initial contact point of the next
            // curve.
            const Vec3& initialContactPoint_S =
                nextCurve->getDataInst(state).X_SP.p();
            // Transform the contact point to ground.
            const Transform& X_GS = nextCurve->calcSurfaceFrameInGround(state);
            return X_GS.shiftFrameStationToBase(initialContactPoint_S);
        };
        // There are no curve segments following given obstacle: the next point
        // is the path's termination point.
        return getTerminationBody()
            .getBodyTransform(state)
            .shiftFrameStationToBase(m_TerminationPoint);
    }
    // Similarly find the first curve segment before the given curve segment.
    const CurveSegment* findPrevActiveCurveSegment(
        const State& state,
        ObstacleIndex ix) const
    {
        for (int i = ix - 1; i >= 0; --i) {
            // Find the active segment before the current.
            if (m_CurveSegments.at(ObstacleIndex(i))
                    .isInContactWithSurface(state)) {
                return &m_CurveSegments.at(ObstacleIndex(i));
            }
        }
        return nullptr;
    }

    // Similarly find the first curve segment after the given curve segment.
    const CurveSegment* findNextActiveCurveSegment(
        const State& state,
        ObstacleIndex ix) const
    {
        // Find the active segment after the current.
        for (int i = ix + 1; i < m_CurveSegments.size(); ++i) {
            if (m_CurveSegments.at(ObstacleIndex(i))
                    .isInContactWithSurface(state)) {
                return &m_CurveSegments.at(ObstacleIndex(i));
            }
        }
        return nullptr;
    }

    // Compute the path error vector (TODO see Scholz2015).
    template <size_t N>
    void calcPathErrorVector(
        const State& state,
        const std::vector<LineSegment>& lines,
        std::array<CoordinateAxis, N> axes,
        Vector& pathError) const;

    // Compute the path error jacobian (TODO see Scholz2015).
    template <size_t N>
    void calcPathErrorJacobian(
        const State& state,
        const std::vector<LineSegment>& lines,
        std::array<CoordinateAxis, N> axes,
        Matrix& J) const;

    //------------------------------------------------------------------------------

    // Reference back to the subsystem.
    CableSubsystem* m_Subsystem = nullptr;
    CableSpanIndex m_Index      = CableSpanIndex::Invalid();

    // TOPOLOGY CACHE (set during realizeTopology())
    CacheEntryIndex indexDataPos = CacheEntryIndex::Invalid();
    CacheEntryIndex indexDataVel = CacheEntryIndex::Invalid();

    MobilizedBodyIndex m_OriginBody = MobilizedBodyIndex::Invalid();
    Vec3 m_OriginPoint{NaN};

    MobilizedBodyIndex m_TerminationBody = MobilizedBodyIndex::Invalid();
    Vec3 m_TerminationPoint{NaN};

    Array_<CurveSegment, ObstacleIndex> m_CurveSegments{};

    Real m_PathAccuracy       = 1e-4;
    int m_SolverMaxIterations = 50; // TODO set to something reasonable.

    // For each curve segment the max allowed radial curvature.
    Real m_MaxCorrectionStepDeg = 10.; // TODO describe

    IntegratorTolerances m_IntegratorTolerances{};

    friend CableSpan;
    friend CableSubsystem;
    friend CableSubsystemTestHelper;
};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

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
    const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);
    const Vec3& x_GB = curve.getMobilizedBody().getBodyOriginLocation(s);

    // Contact point moment arms in ground.
    const Vec3 r_P = dataPos.X_GP.p() - x_GB;
    const Vec3 r_Q = dataPos.X_GQ.p() - x_GB;

    // Tangent directions at contact points in ground.
    const UnitVec3& t_P = getTangent(dataPos.X_GP);
    const UnitVec3& t_Q = getTangent(dataPos.X_GQ);

    unitForce_G[0] = r_Q % t_Q - r_P % t_P;
    unitForce_G[1] = t_Q - t_P;
}

void calcUnitForceAtCableOrigin(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    // Origin contact point moment arm in ground.
    const CableSpanData::Pos dataPos = cable.getPosInfo(s);
    const Vec3& arm_G =
        dataPos.originPoint_G - cable.getOriginBody().getBodyOriginLocation(s);

    unitForce_G[0] = -arm_G % dataPos.originTangent_G;
    unitForce_G[1] = -Vec3(dataPos.originTangent_G);
}

void calcUnitForceAtCableTermination(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CableSpanData::Pos dataPos = cable.getPosInfo(s);
    const Vec3& arm_G                = dataPos.terminationPoint_G -
                        cable.getTerminationBody().getBodyOriginLocation(s);

    unitForce_G[0] = -arm_G % dataPos.terminationTangent_G;
    unitForce_G[1] = -Vec3(dataPos.terminationTangent_G);
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

    const bool isInContactWithSurface    = curve.isInContactWithSurface(s);
    const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);

    // Draw the surface (TODO since we do not own it, should it be done here at
    // all?).
    {
        DecorativeGeometry geo = curve.getDecoration(); // TODO clone it?
        const Transform& X_GS  = dataPos.X_GS;
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
        auto drawLine = [&](Real length, Vec3 nextPoint) {
            if (length != 0.) {
                decorations.push_back(DecorativeLine(prevPoint, nextPoint)
                                          .setColor(Purple)
                                          .setLineThickness(3));
            }
            prevPoint = nextPoint;
        };

        // TODO for debugging: Choose to show the original intergrator's path
        // TODO this is actually broken
        // points or the resampled path points.
        constexpr bool doResampling = false;
        if (doResampling) {
            curve.calcPathPoints(s, drawLine);
        } else {
            // Use points from GeodesicIntegrator directly.
            for (const ContactGeometry::ImplicitGeodesicState& q :
                 curve.getDataInst(s).geodesicIntegratorStates) {
                drawLine(
                    q.arcLength,
                    dataPos.X_GS.shiftFrameStationToBase(q.point));
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

        DrawFrenetFrame(dataPos.X_GP);
        DrawFrenetFrame(dataPos.X_GQ);
    }

    // Draw the initial contact point hint using a yellow line.
    // TODO this is for debugging should be removed.
    {
        constexpr int c_HintLineThickness = 2;

        const Transform& X_GS = dataPos.X_GS;

        decorations.push_back(DecorativeLine(
                                  dataPos.X_GS.p(),
                                  dataPos.X_GS.shiftFrameStationToBase(
                                      curve.getContactPointHint()))
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

    const CableSpanData::Pos& dataPos = getPosInfo(s);
    Vec3 prevPoint                    = dataPos.originPoint_G;

    for (const CurveSegment& curve : m_CurveSegments) {
        if (curve.isInContactWithSurface(s)) {
            const CurveSegmentData::Pos cppe = curve.getDataPos(s);

            const Vec3 nextPoint = cppe.X_GP.p();
            decorations.push_back(DecorativeLine(prevPoint, nextPoint)
                                      .setColor(Purple)
                                      .setLineThickness(c_LineThickness));
            prevPoint = cppe.X_GQ.p();
        }

        calcCurveDecorativeGeometryAndAppend(curve, s, decorations);
    }

    // TODO choose colors.
    decorations.push_back(DecorativeLine(prevPoint, dataPos.terminationPoint_G)
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
            const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);

            // Compute the line segment from the previous point to the
            // first curve contact point.
            lines.emplace_back(prevPathPoint, dataPos.getFirstContactPoint());

            // The next line segment will start at the curve's final contact
            // point.
            prevPathPoint = dataPos.getFinalContactPoint();
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
            totalCableLength += curve.getDataInst(s).length;
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
    int row    = -1;

    // Reset path error vector to zero.
    pathError *= 0;

    forEachActiveCurveSegment(
        *this,
        s,
        // Compute the path error for each curved segment.
        [&](const CurveSegment& curve) {
            const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);
            // Compute path error at first contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++row) =
                    calcPathError(lines.at(lineIx), dataPos.X_GP.R(), axis);
            }
            ++lineIx;
            // Compute path error at final contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++row) =
                    calcPathError(lines.at(lineIx), dataPos.X_GQ.R(), axis);
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
    const Transform& X_GP = curve.getDataPos(state).X_GP;

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
    const Transform& X_GQ = curve.getDataPos(s).X_GQ;

    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);
    const Real a            = dataInst.jacobi_Q[0];
    const Real r            = dataInst.jacobi_Q[1];

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
    const Transform& X_GP = curve.getDataPos(s).X_GP;

    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);

    const Real tau = dataInst.torsion_P;
    const Real kt  = dataInst.curvatures_P[0];
    const Real kb  = dataInst.curvatures_P[1];

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
    const Transform& X_GQ = curve.getDataPos(s).X_GQ;

    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);

    const Real tau = dataInst.torsion_Q;
    const Real kt  = dataInst.curvatures_Q[0];
    const Real kb  = dataInst.curvatures_Q[1];

    const Real a = dataInst.jacobi_Q[0];
    const Real r = dataInst.jacobi_Q[1];

    const Real aDot = dataInst.jacobiDot_Q[0];
    const Real rDot = dataInst.jacobiDot_Q[1];

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

        const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);
        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_P = dataPos.X_GP.R().getAxisUnitVec(axis);

            AddBlock(calcJacobianOfPathErrorAtP(curve, s, l_P, a_P));

            if (prev) {
                AddBlock(calcJacobianOfNextPathError(*prev, s, l_P, a_P), -Nq);
            }
            ++row;
        }

        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_Q = dataPos.X_GQ.R().getAxisUnitVec(axis);

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
// jacobian in the MatrixWorkspace. The result is a vector of Corrections for
// each curve.
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

    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);

    // Clamp tangential displacement at the initial contact point.
    {
        const Real dxEst = c[0];
        const Real k     = dataInst.curvatures_P[0];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp binormal displacement at the initial contact point.
    {
        const Real dxEst = c[1];
        const Real k     = dataInst.curvatures_P[1];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp tangential displacement at the final contact point.
    {
        const Real dxEst = std::abs(c[0]) + std::abs(c[3]);
        const Real k     = dataInst.curvatures_Q[0];
        UpdateMaxStepSize(dxEst, k);
    }

    // Clamp binormal displacement at the final contact point.
    {
        const Real a     = dataInst.jacobi_Q[0];
        const Real r     = dataInst.jacobi_Q[1];
        const Real dxEst = std::abs(c[1] * a) + std::abs(c[2] * r);
        const Real k     = dataInst.curvatures_Q[1];
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
    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);
    const FrenetFrame& X_SP    = dataInst.X_SP;

    const Real tau   = dataInst.torsion_P;
    const Real kappa = dataInst.curvatures_P[0];

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
        std::max(dataInst.length + dl, 0.); // Clamp length to be nonnegative.

    // Shoot the new geodesic.
    curve.calcLocalGeodesic(
        s,
        correctedPoint_S,
        correctedTangent_S,
        correctedLength,
        dataInst.integratorInitialStepSize,
        tols);
}

} // namespace

//==============================================================================
//                      CABLE SPAN POSITION LEVEL CACHE
//==============================================================================

void CableSpan::Impl::calcPosInfo(const State& s, CableSpanData::Pos& dataPos)
    const
{
    // Path origin and termination point.
    const Vec3 x_O =
        getOriginBody().getBodyTransform(s).shiftFrameStationToBase(
            m_OriginPoint);
    const Vec3 x_I =
        getTerminationBody().getBodyTransform(s).shiftFrameStationToBase(
            m_TerminationPoint);

    dataPos.originPoint_G      = x_O;
    dataPos.terminationPoint_G = x_I;

    // Axes considered when computing the path error.
    const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

    dataPos.loopIter = 0;
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
            dataPos.cableLength          = (x_I - x_O).norm();
            dataPos.terminationTangent_G = dataPos.originTangent_G =
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
        dataPos.pathError = data.pathError.normInf();

        // Stop iterating if max path error is small, or max iterations has been
        // reached.
        if (dataPos.pathError < m_PathAccuracy ||
            dataPos.loopIter >= m_SolverMaxIterations) {
            // Update cache entry and stop solver.
            dataPos.cableLength =
                calcTotalCableLength(*this, s, data.lineSegments);
            dataPos.originTangent_G      = data.lineSegments.front().direction;
            dataPos.terminationTangent_G = data.lineSegments.back().direction;
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
        calcPathCorrections(data, dataPos.pathError);

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

        ++dataPos.loopIter;
    }
}

//==============================================================================
//                      CABLE SPAN VELOCITY LEVEL CACHE
//==============================================================================

void CableSpan::Impl::calcVelInfo(const State& s, CableSpanData::Vel& velInfo)
    const
{
    const CableSpanData::Pos& dataPos = getPosInfo(s);

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
        const MobilizedBody& mobod                = curve.getMobilizedBody();
        const CurveSegmentData::Pos& curveDataPos = curve.getDataPos(s);

        const UnitVec3& e_G = getTangent(curveDataPos.X_GP);

        const Vec3 v_GP =
            CalcPointVelocityInGround(mobod, curveDataPos.X_GP.p());

        lengthDot += dot(e_G, v_GP - v_GQ);

        x_GQ = curveDataPos.X_GQ.p();
        v_GQ = CalcPointVelocityInGround(mobod, x_GQ);
    });

    const Vec3 v_GP =
        getTerminationBody().findStationVelocityInGround(s, m_TerminationPoint);

    const UnitVec3 e_G(dataPos.terminationPoint_G - x_GQ);

    lengthDot += dot(e_G, v_GP - v_GQ);
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
    CurveSegmentData::Instance& dataInst)
{
    dataInst.geodesicIntegratorStates.clear();
    dataInst.integratorInitialStepSize = initIntegratorStepSize;

    const ContactGeometry& geometry = curve.getContactGeometry();

    // TODO implement shooting Analytic and Parametric geodesics.
    geometry.shootGeodesicInDirectionImplicitly(
        point_S,
        tangent_S,
        length,
        initIntegratorStepSize,
        tols.intergatorAccuracy,
        tols.constraintProjectionTolerance,
        tols.constraintProjectionMaxIterations,
        [&](const ContactGeometry::ImplicitGeodesicState& q) {
            // Log the integrator knot points.
            dataInst.geodesicIntegratorStates.push_back(q);

            dataInst.jacobi_Q    = {q.jacobiTrans, q.jacobiRot};
            dataInst.jacobiDot_Q = {q.jacobiTransDot, q.jacobiRotDot};
        });

    // Update the initial integrator step size.
    {
        // We want to know if we should attempt a larger initIntegratorStepSize
        // next time, or a smaller one. If the initIntegratorStepSize was
        // rejected, it makes little sense to start with the exact same step
        // size next time, only to find it rejected again.

        // If the integrator took a single step or none: The final arc length is
        // shorter than the initial step size attempted. The step was atleast
        // not rejected, but we do not know if we can make it larger. So we do
        // not update it.
        const int numberOfIntegrationSteps =
            dataInst.geodesicIntegratorStates.size() - 1;
        if (numberOfIntegrationSteps <= 1) {
            // Try the same step next time:
            dataInst.integratorInitialStepSize = initIntegratorStepSize;
        }

        // If the integrator took two or more steps, we might want to update the
        // next step size to try,
        if (numberOfIntegrationSteps >= 2) {
            // The acual initial step size taken by the integrator cannot have
            // been larger than initIntegratorStepSize, but it might be smaller
            // if it was rejected. If rejected, we will reduce the init step
            // size for next time.
            const Real actualFirstStepSize =
                dataInst.geodesicIntegratorStates.at(1).arcLength -
                dataInst.geodesicIntegratorStates.front().arcLength;
            const bool initStepWasRejected =
                actualFirstStepSize < initIntegratorStepSize;
            if (initStepWasRejected) {
                // Reduce the init step size for next time, to avoid rejecting
                // it again.
                dataInst.integratorInitialStepSize = actualFirstStepSize;
            } else {
                // If the initIntegratorStepSize was accepted we might want to
                // try a larger step next time. The actualFirstStepSize cannot
                // be larger than initIntegratorStepSize: we must take a look
                // at the second step size.
                const Real actualSecondStepSize =
                    dataInst.geodesicIntegratorStates.at(2).arcLength -
                    dataInst.geodesicIntegratorStates.at(1).arcLength;
                // We already established that the initIntegratorStepSize was
                // accepted, now we are only looking to increase the step size.
                dataInst.integratorInitialStepSize =
                    std::max(initIntegratorStepSize, actualSecondStepSize);
            }
        }
    }

    dataInst.X_SP = calcFrenetFrameFromGeodesicState(
        geometry,
        dataInst.geodesicIntegratorStates.front());
    dataInst.X_SQ = calcFrenetFrameFromGeodesicState(
        geometry,
        dataInst.geodesicIntegratorStates.back());

    dataInst.curvatures_P = {
        calcSurfaceCurvature(geometry, dataInst.X_SP, TangentAxis),
        calcSurfaceCurvature(geometry, dataInst.X_SP, BinormalAxis),
    };

    dataInst.curvatures_Q = {
        calcSurfaceCurvature(geometry, dataInst.X_SQ, TangentAxis),
        calcSurfaceCurvature(geometry, dataInst.X_SQ, BinormalAxis),
    };

    dataInst.torsion_P = geometry.calcSurfaceTorsionInDirection(
        dataInst.X_SP.p(),
        getTangent(dataInst.X_SP));
    dataInst.torsion_Q = geometry.calcSurfaceTorsionInDirection(
        dataInst.X_SQ.p(),
        getTangent(dataInst.X_SQ));

    dataInst.length = length;
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
    if (curve.getDataInst(s).status != WrappingStatus::InitialGuess) {
        return;
    }
    // The first integrator stepsize to attempt is a bit arbitrary: Taking it
    // too large will simply cause it to be rejected. Taking it too small will
    // cause it to rapidly increase. We take it too small, to be on the safe
    // side.
    constexpr Real c_initIntegratorStepSize = 1e-10;

    // Shoot a zero length geodesic at the contact point, with the tangent
    // directed from prevPoint to nextPoint.
    curve.calcLocalGeodesic(
        s,
        curve.getContactPointHint(), // = location
        nextPoint_S - prevPoint_S,   // = direction
        0.,                          // = arc length
        c_initIntegratorStepSize,
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
    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);
    if (dataInst.status != WrappingStatus::LiftedFromSurface) {
        return;
    }

    // Use the cached tracking point as the initial guess.
    Vec3 pointOnLineNearSurface_S = dataInst.trackingPointOnLine_S;

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
            dataInst.integratorInitialStepSize,
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
    const CurveSegmentData::Instance& dataInst = curve.getDataInst(s);
    if (dataInst.status != WrappingStatus::InContactWithSurface) {
        return;
    }

    // The curve length must have shrunk completely before lifting off.
    if (dataInst.length > 0.) {
        return;
    }

    // For a zero-length curve, trigger liftoff when the prev and next points
    // lie above the surface plane.
    if (dot(prevPoint_S - dataInst.X_SP.p(),
            dataInst.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0. ||
        dot(nextPoint_S - dataInst.X_SP.p(),
            dataInst.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0.) {
        // No liftoff.
        return;
    }

    // Liftoff detected, initialize the tracking point from the last contact
    // point.
    curve.liftCurveFromSurface(s, dataInst.X_SP.p());
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
    CurveSegmentData::Instance& dataInst = updDataInst(s);

    dataInst.status = WrappingStatus::InContactWithSurface;

    shootNewGeodesic(
        *this,
        point_S,
        tangent_S,
        length,
        stepSizeHint,
        tols,
        dataInst);

    getSubsystem().markDiscreteVarUpdateValueRealized(s, m_InstanceIx);

    invalidatePosEntry(s);
}

void CurveSegment::liftCurveFromSurface(const State& s, Vec3 trackingPoint_S)
    const
{
    CurveSegmentData::Instance& dataInst = updDataInst(s);

    dataInst.status                = WrappingStatus::LiftedFromSurface;
    dataInst.trackingPointOnLine_S = trackingPoint_S;

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

    if (getDataInst(s).status == WrappingStatus::Disabled) {
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
    const CurveSegmentData::Instance& dataInst = getDataInst(s);

    // Start updating the position level cache.
    CurveSegmentData::Pos& dataPos = updPosInfo(s);

    // Transform geodesic in local surface coordinates to ground.
    {
        // Store the local geodesic in ground frame.
        dataPos.X_GS = X_GS;

        // Store the local geodesic in ground frame.
        dataPos.X_GP = X_GS.compose(dataInst.X_SP);
        dataPos.X_GQ = X_GS.compose(dataInst.X_SQ);
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
    const ContactGeometry::ImplicitGeodesicState& a,
    const ContactGeometry::ImplicitGeodesicState& b,
    Real l)
{
    return calcHermiteInterpolation(
        a.arcLength,
        a.point,
        a.tangent,
        b.arcLength,
        b.point,
        b.tangent,
        l);
}

// Compute the tangent in between two geodesic samples using Hermite
// interpolation.
Vec3 calcInterpolatedTangent(
    const ContactGeometry& geometry,
    const ContactGeometry::ImplicitGeodesicState& a,
    const ContactGeometry::ImplicitGeodesicState& b,
    Real l)
{
    // Helper for calculating the derivative of the tangent.
    auto CalcTangentDot =
        [&](const ContactGeometry::ImplicitGeodesicState& q) -> Vec3 {
        const FrenetFrame& X_S = calcFrenetFrameFromGeodesicState(geometry, q);
        const Real k = calcSurfaceCurvature(geometry, X_S, TangentAxis);
        return k * getNormal(X_S);
    };

    return calcHermiteInterpolation(
        a.arcLength,
        a.tangent,
        CalcTangentDot(a),
        b.arcLength,
        b.tangent,
        CalcTangentDot(b),
        l);
}

using InterpolatorFn = std::function<void(
    const ContactGeometry::ImplicitGeodesicState& a,
    const ContactGeometry::ImplicitGeodesicState& b,
    Real l)>;

// Resamples the geodesic samples at equal length intervals.
// @param geodesic contains the original geodesic samples.
// @param nSamples the number of samples to obtain.
// @param interpolator function for performing the interpolation between two
// selected samples at a given length.
// @return the number of samples used for the interpolation.
int resampleGeodesic(
    const std::vector<ContactGeometry::ImplicitGeodesicState> geodesic,
    int nSamples, // TODO use length increment?
    InterpolatorFn& interpolator)
{
    // Some sanity checks.
    if (geodesic.empty()) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Resampling of geodesic failed: Provided geodesic is empty.");
    }
    if (geodesic.front().arcLength != 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: First frame "
                                 "must be at arcLength = zero");
    }
    if (geodesic.back().arcLength < 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: Last frame "
                                 "must be at arcLength > zero");
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
    const Real dl = geodesic.back().arcLength / static_cast<Real>(nSamples - 1);

    // Compute the interpolated points from the geodesic.
    auto itGeodesic = geodesic.begin();
    // We can skip the first and last samples, because these are pushed
    // manually before and after this loop respectively (we start at i=1
    // and stop at i < nSamples-1).
    for (int i = 0; i < nSamples - 1; ++i) {

        // Length at the current interpolation point.
        const Real arcLength = dl * static_cast<Real>(i);

        // Find the two samples (lhs, rhs) of the geodesic such that the
        // arcLength of the interpolation point lies between them.
        // i.e. find: lhs.arcLength <= arcLength < rhs.arcLength
        while (true) {
            // Sanity check: We should stay within range.
            if ((itGeodesic + 1) == geodesic.end()) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Attempted to read out of "
                    "array range");
            }

            // The candidate samples to use for interpolation.
            const ContactGeometry::ImplicitGeodesicState& lhs = *itGeodesic;
            const ContactGeometry::ImplicitGeodesicState& rhs =
                *(itGeodesic + 1);

            // Sanity check: Samples are assumed to be monotonically increasing
            // in arcLength.
            if (lhs.arcLength > rhs.arcLength) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Samples are not "
                    "monotonically increasing in arcLength.");
            }

            // Check that the interpolation point lies between these samples:
            // lhs.arcLength <= arcLength < rhs.arcLength
            if (arcLength >= rhs.arcLength) {
                // Try the next two samples.
                ++itGeodesic;
                continue;
            }

            // Write interpolated point to the output buffer.
            interpolator(lhs, rhs, arcLength);

            break;
        }
    }

    // Capture the last point of the geodesic.
    interpolator(geodesic.front(), geodesic.back(), geodesic.back().arcLength);

    return nSamples;
}
} // namespace

int CurveSegment::calcPathPointsAndTangents(
    const State& state,
    const std::function<void(Real length, Vec3 point_G, UnitVec3 tangent_G)>&
        sink,
    int nSamples) const
{
    if (!isInContactWithSurface(state)) {
        return 0;
    }

    const CurveSegmentData::Instance& dataInst = getDataInst(state);
    const CurveSegmentData::Pos& dataPos = getDataPos(state);

    // If the curve has zero length, return a single point, nothing to
    // interpolate.
    if (dataInst.length == 0.) {
        sink(0., dataPos.X_GP.p(), getTangent(dataPos.X_GP));
        return 1;
    }

    // Requesting two inerpolation points means we write the first and last
    // contact point. Nothing to interpolate.
    if (nSamples == 2) {
        sink(0., dataPos.X_GP.p(), getTangent(dataPos.X_GP));
        sink(dataInst.length, dataPos.X_GQ.p(), getTangent(dataPos.X_GQ));
        return 2;
    }

    // Resample the points from the integrator by interpolating at equal length
    // intervals.
    InterpolatorFn Interpolator =
        [&](const ContactGeometry::ImplicitGeodesicState& a,
            const ContactGeometry::ImplicitGeodesicState& b,
            Real l) {
            const Transform& X_GS = dataPos.X_GS;
            const Vec3 interpolatedPoint_G =
                X_GS.shiftFrameStationToBase(calcInterpolatedPoint(a, b, l));
            const UnitVec3 interpolatedTangent_G =
                UnitVec3(X_GS.xformFrameVecToBase(
                    calcInterpolatedTangent(getContactGeometry(), a, b, l)));

            sink(l, interpolatedPoint_G, interpolatedTangent_G);
        };
    return resampleGeodesic(
        dataInst.geodesicIntegratorStates,
        nSamples,
        Interpolator);
}

int CurveSegment::calcPathPoints(
    const State& state,
    const std::function<void(Real length, Vec3 point_G)>& sink,
    int nSamples) const
{
    if (!isInContactWithSurface(state)) {
        return 0;
    }

    const CurveSegmentData::Instance& dataInst = getDataInst(state);
    const CurveSegmentData::Pos& dataPos = getDataPos(state);

    // If the curve has zero length, return a single point, nothing to
    // interpolate.
    if (dataInst.length == 0.) {
        sink(0., dataPos.X_GP.p());
        return 1;
    }

    // Requesting two inerpolation points means we write the first and last
    // contact point. Nothing to interpolate.
    if (nSamples == 2) {
        sink(0., dataPos.X_GP.p());
        sink(dataInst.length, dataPos.X_GQ.p());
        return 2;
    }

    // Resample the points from the integrator by interpolating at equal
    // intervals.

    // Define the required function that interpolates between two samples, and
    // logs the interpolated sample.
    InterpolatorFn Interpolator =
        [&](const ContactGeometry::ImplicitGeodesicState& a,
            const ContactGeometry::ImplicitGeodesicState& b,
            Real l) {
            const Vec3 interpolatedPoint_S = calcInterpolatedPoint(a, b, l);
            const Vec3 interpolatedPoint_G =
                dataPos.X_GS.shiftFrameStationToBase(interpolatedPoint_S);
            sink(l, interpolatedPoint_G);
        };
    return resampleGeodesic(
        dataInst.geodesicIntegratorStates,
        nSamples,
        Interpolator);
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

        // Count the number of active curve segments.
        const int nActive = cable.countActive(s);

        if (nActive == 0) {
            os << "Cable has no active segments: Skipping perturbation test\n";
            continue;
        }

        constexpr int DOF = 4;
        for (int i = 0; i < (DOF * 2 * nActive); ++i) {
            // We do not want to mess with the actual state, so we make a copy.
            const State sCopy = s;
            system.realize(sCopy, Stage::Position);

            // Trigger realizing position level cache, resetting the
            // configuration.
            const CableSpanData::Pos dataPos = cable.getPosInfo(sCopy);

            SimTK_ASSERT(
                cable.countActive(sCopy) == nActive,
                "Unexpected change in number of wrapping segments during "
                "perturbation test");

            MatrixWorkspace& data =
                subsystem.getImpl().updSolverData(sCopy).updOrInsert(nActive);

            // Define the perturbation we will use for testing the jacobian.
            perturbation *= -1.;
            data.pathCorrection.setTo(0.);
            data.pathCorrection.set(i / 2, perturbation);

            calcLineSegments(
                cable,
                sCopy,
                dataPos.originPoint_G,
                dataPos.terminationPoint_G,
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
                prevWrappingStatus.push_back(
                    curve.getDataInst(sCopy).status);
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
                dataPos.originPoint_G,
                dataPos.terminationPoint_G,
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
                    WrappingStatusChanged |=
                        prevWrappingStatus.at(++ix) !=
                        curve.getDataInst(sCopy).status;
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

//==============================================================================
//                          CableSubsystem::Impl
//==============================================================================

int CableSubsystem::Impl::realizeSubsystemTopologyImpl(State& state) const
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

const CableSpan& CableSubsystem::Impl::getCable(CableSpanIndex index) const
{
    SimTK_ERRCHK1_ALWAYS(
        &cables[index].getImpl().getSubsystem().getImpl() == this,
        "CableSubsystem::getCable",
        "Cable %d is not owned by this subsystem",
        index);
    SimTK_ERRCHK1_ALWAYS(
        cables[index].getImpl().getIndex() == index,
        "CableSubsystem::getCable",
        "Cable %d has an invalid index",
        index);
    return cables[index];
}

CableSpan& CableSubsystem::Impl::updCable(CableSpanIndex index)
{
    SimTK_ERRCHK1_ALWAYS(
        &cables[index].getImpl().getSubsystem().getImpl() == this,
        "CableSubsystem::updCable",
        "Cable %d is not owned by this subsystem",
        index);
    SimTK_ERRCHK1_ALWAYS(
        cables[index].getImpl().getIndex() == index,
        "CableSubsystem::getCable",
        "Cable %d has an invalid index",
        index);
    return cables[index];
}

int CableSubsystem::Impl::calcDecorativeGeometryAndAppendImpl(
    const State& state,
    Stage stage,
    Array_<DecorativeGeometry>& decorations) const
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

//==============================================================================
//                          Cable Subsystem
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

const MultibodySystem& CableSubsystem::getMultibodySystem() const
{
    return getImpl().getMultibodySystem();
}

//==============================================================================
//                                Cable Span
//==============================================================================

// Allocate a new default-constructed CableSpan::Impl.
CableSpan::CableSpan() : m_Impl(new Impl())
{}

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
    CableSpanIndex ix = subsystem.updImpl().adoptCable(*this);
    updImpl().setSubsystem(subsystem, ix);
}

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
void CableSpan::setObstacleXformSurfaceToBody(
    ObstacleIndex ix,
    const Transform& X_BS)
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
    return getImpl().getCurveSegment(ix).getDataInst(state).status;
}

Real CableSpan::getCurveSegmentLength(const State& state, ObstacleIndex ix)
    const
{
    getImpl().realizePosition(state);
    return getImpl().getCurveSegment(ix).getDataInst(state).length;
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
    const std::function<void(Real length, Vec3 point, UnitVec3 tangent)>& sink)
    const
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
    return getImpl()
        .getCurveSegment(ix)
        .getDataInst(state)
        .geodesicIntegratorStates.size();
}

Real CableSpan::getCurveSegmentInitialIntegratorStepSize(
    const State& state,
    ObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl()
        .getCurveSegment(ix)
        .getDataInst(state)
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
    for (ObstacleIndex ix(0); ix < getImpl().getNumCurveSegments(); ++ix) {
        getImpl().getCurveSegment(ix).storeCurrentPath(state);
    }
}
