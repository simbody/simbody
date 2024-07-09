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

#include "SimTKcommon/internal/State.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

using namespace SimTK;

// For convenience.
using ObstacleIndex = CableSpanObstacleIndex;

//==============================================================================
//                            Data Structures
//==============================================================================
/* The following section contains the data structures used to compute the cable
path. */
namespace
{

//------------------------------------------------------------------------------
//  Geodesic Degrees of Freedom
//------------------------------------------------------------------------------
/* Number of degrees of freedom of a general geodesic.
There are four: see NaturalGeodesicCorrections below. */
constexpr int c_GeodesicDOF = 4;

//------------------------------------------------------------------------------
//  Natural Geodesic Corrections
//------------------------------------------------------------------------------
/* Natural geodesic corrections represent the coordinates along which we can
change a geodesic. The elements, in order, are (see Scholz2015):
- Tangential correction of initial contact point.
- Binormal correction of initial contact point.
- Directional correction of tangent at the initial contact point.
- Lengthening correction of geodesic. */
using NaturalGeodesicCorrection = Vec<c_GeodesicDOF>;

//------------------------------------------------------------------------------
//  Struct Line Segment
//------------------------------------------------------------------------------
/* Representation of a cable segment that does not lie on a surface: A straight
line. */
struct LineSegment final {
    LineSegment() = default;

    // Construct a new LineSegment connecting pointA and pointB.
    LineSegment(const Vec3& pointA, const Vec3& pointB) :
        length((pointB - pointA).norm()), direction((pointB - pointA) / length)
    {}

    Real length = NaN;
    UnitVec3 direction{NaN, NaN, NaN};
};

//------------------------------------------------------------------------------
//  Wrapping Status
//------------------------------------------------------------------------------
/* This enum gives the wrapping status of an obstacle in a CableSpan's path. */
enum class ObstacleWrappingStatus {
    // InitialGuess: The curve segment must be initialized from the initial
    // contact point parameter. This will be used as the starting point for
    // computing the optimal cable path.
    InitialGuess,
    // InContactWithSurface: The cable comes into contact with the obstacle.
    InContactWithSurface,
    // LiftedFromSurface: The cable does not come into contact with the
    // obstacle.
    LiftedFromSurface,
    // Disabled: The obstacle was manually disabled, preventing any interaction
    // with the cable. TODO no longer used?
    Disabled,
};

//------------------------------------------------------------------------------
//  Frenet Frame
//------------------------------------------------------------------------------
/* The frenet frame plays an important role in defining the state of a geodesic.
In the end it is simply a Transform, which is why we define it as an alias.

The position of the frenet frame is defined to lie on the geodesic, but for
the orientation of the frame different conventions are used.
Here we define it as (see Scholz2015):
- X axis: tangent to geodesic
- Y axis: normal to surface
- Z axis: binormal to geodesic (tangent cross normal) */
using FrenetFrame                        = Transform;
static const CoordinateAxis TangentAxis  = XAxis;
static const CoordinateAxis NormalAxis   = YAxis;
static const CoordinateAxis BinormalAxis = ZAxis;

//------------------------------------------------------------------------------
//  Struct MatrixWorkspace
//------------------------------------------------------------------------------
/* This is a helper struct that is used by a CableSpan to compute the
Stage::Position level data, i.e. the spanned path.
Computing the path involves a Newton type iteration, for which several
matrices need to be computed. This struct provides the required matrices with
appropriate dimensions. After computing the path this data is no longer needed
by the CableSpan. */
struct MatrixWorkspace {
    // Number of path error constraints per curve segment.
    static constexpr int c_NumPathErrorConstraints = 4;
    // Total number of constraints per curve segment.
    static constexpr int c_NumConstraints = c_NumPathErrorConstraints + 1;

    // Given the number of CurveSegments that are in contact with their
    // respective obstacle's surface, contruct a MatrixWorkspace of correct
    // dimensions.
    explicit MatrixWorkspace(int problemSize) : nObstaclesInContact(problemSize)
    {
        static constexpr int Q = c_GeodesicDOF;
        static constexpr int C = c_NumConstraints;
        const int n            = problemSize;

        lineSegments.resize(n + 1);
        pathErrorJacobian = Matrix(C * n, Q * n, 0.);
        pathCorrection    = Vector(Q * n, 0.);
        pathError         = Vector(C * n, 0.);
    }

    // Return the NaturalGeodesicCorrection for the curve segment at the
    // "active" curve index, where "active" means counting those that are in
    // contact with the obstacle.
    NaturalGeodesicCorrection getCurveCorrection(int activeCurveIx) const
    {
        SimTK_ASSERT(
            activeCurveIx >= nObstaclesInContact,
            "Index of curve (counting active only) exceeds number of segments in contact with obstacles");
        SimTK_ASSERT(
            pathCorrection.size() == nObstaclesInContact * c_GeodesicDOF,
            "Invalid size of pathCorrection vector.");
        int eltIx = activeCurveIx * c_GeodesicDOF;
        return {
            pathCorrection[eltIx],
            pathCorrection[eltIx + 1],
            pathCorrection[eltIx + 2],
            pathCorrection[eltIx + 3]};
    }

    // All straight line segments in the current path.
    std::vector<LineSegment> lineSegments;

    // TODO describe
    Vector pathError;
    Matrix pathErrorJacobian;
    Vector pathCorrection;
    FactorQTZ inverse;

    Real maxPathError       = NaN;
    int nObstaclesInContact = -1;
};

//------------------------------------------------------------------------------
//  CableSpanData
//------------------------------------------------------------------------------
/* CableSpanData is a data structure used by class CableSpan::Impl to store
relevant quantities in the State's cache.

The struct has three members Instance, Pos and Vel that correspond to the
stages (Stage::Instance, Stage::Position and Stage::Velocity) at which they can
be computed.

Member variables with a "_G" suffix are expressed in the Ground frame. */
struct CableSpanData {
    /* Cache entry for holding MatrixWorkspaces of different dimensions.
    This acts as a scratchpad when computing Stage::Position level data, during
    which the path is essentially computed. During that computation the
    dimension of the matrices changes as the cable touches down- or lifts off
    from obstacle surfaces. This means that MatrixWorkspaces of different
    dimensions are used when compting the optimal path. After having completed
    the Stage::Position level data, the data in the MatrixWorkspaces is no
    longer used. */
    class Instance {
    public:
        Instance() = default;

        // Get mutable access to a MatrixWorkspace of appropriate dimension.
        // Constructs requested MatrixWorkspace if not done previously.
        MatrixWorkspace& updOrInsert(int nActive)
        {
            SimTK_ASSERT1(
                nActive >= 0,
                "CableSpanData::updOrInsert()"
                "Number of obstacles in contact must be nonnegative (got %d)",
                nActive);

            // Construct all MatrixWorkspaces for up to and including the
            // requested dimension.
            for (size_t i = matrixWorkspaces.size(); i <= nActive; ++i) {
                matrixWorkspaces.emplace_back(i);
            }

            // Return MatrixWorkspace of requested dimension.
            return matrixWorkspaces.at(nActive);
        }

    private:
        // MatrixWorkspaces of suitable dimensions for solving the cable path
        // given a number of obstacles in contact with the cable. The ith
        // element in the vector is used for solving a path with i-active
        // obstacles.
        std::vector<MatrixWorkspace> matrixWorkspaces;
    };
    struct Pos final {
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
        // The path smoothness, computed as the maximum angular misalignment at
        // the points where the straight- and curved-line segments meet,
        // measured in radians.
        Real smoothness = NaN;
        // Number of iterations the solver used to compute the cable's path.
        int loopIter = -1;
    };
    struct Vel final {
        // Time derivative of the total cable length.
        Real lengthDot = NaN;
    };
};

//------------------------------------------------------------------------------
//  CurveSegmentData
//------------------------------------------------------------------------------
/* CurveSegmentData is a data structure used by class CurveSegment to store
relevent quantities in the State's cache.

The struct has two member structs Instance and Pos that correspond to the stages
(Stage::Instance and Stage::Position) at which they can be computed.

Following Scholz2015 the following subscripts are used:

Member variables with a "_P" suffix indicate the start of the geodesic.
Member variables with a "_Q" suffix indicate the end of the geodesic.
Member variables with a "_S" suffix are expressed in the Surface frame.
Member variables with a "_G" suffix are expressed in the Ground frame.

For example: point_SP would indicate the initial contact point represented and
measured from the surface's origin frame. */
struct CurveSegmentData {
    struct Instance final {
        // Frenet frame at the initial contact point w.r.t. the surface frame.
        FrenetFrame X_SP{Rotation().setRotationToNaN(), Vec3(NaN)};
        // Frenet frame at the final contact point w.r.t. the surface frame.
        FrenetFrame X_SQ{Rotation().setRotationToNaN(), Vec3(NaN)};
        // Length of this curve segment.
        Real arcLength = NaN;
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
        std::vector<ContactGeometry::GeodesicKnotPoint>
            geodesicIntegratorStates;
        // The initial integrator stepsize to try next time when shooting a
        // geodesic. This step size estimate will improve each time after
        // shooting a new geodesic.
        // The default value is a bit arbitrary: Taking it too large will
        // simply cause it to be rejected. Taking it too small will cause it to
        // rapidly increase. We take it too small, to be on the safe side.
        Real integratorInitialStepSize = 1e-10;
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
        ObstacleWrappingStatus wrappingStatus =
            ObstacleWrappingStatus::InitialGuess;
    };
    struct Pos final {
        // Position and orientation of contact geometry w.r.t. ground.
        Transform X_GS{};
        // Frenet frame at the initial contact point w.r.t. ground.
        FrenetFrame X_GP{};
        // Frenet frame at the final contact point w.r.t. ground.
        FrenetFrame X_GQ{};
    };
};

} // namespace

//==============================================================================
//                          Configuration Parameters
//==============================================================================
/* The following section contains structs that collect relevant configuration
parameters. */
namespace
{

//------------------------------------------------------------------------------
//  Integrator Tolerances
//------------------------------------------------------------------------------
/* These are all the parameters needed to configure the GeodesicIntegrator
when computing a new geodesic over the obstacle's surface.

TODO These need reasonable default values. */
struct IntegratorTolerances {
    // The accuracy used to control the stepsize of the variable step geodesic
    // integrator.
    Real geodesicIntegratorAccuracy = 1e-6;
    // TODO the surface projection tolerance should be linked to the geodesic
    // integrator accuracy.
    Real constraintProjectionTolerance = 1e-6;
    // TODO should be used during numerical integration as well, but not
    // connected yet.
    int constraintProjectionMaxIterations = 50;
};

//------------------------------------------------------------------------------
//  Cable Span Parameters
//------------------------------------------------------------------------------
/* These are all the parameters for configuring the solver steps in the
CableSpan. */
struct CableSpanParameters final : IntegratorTolerances {
    // The path's smoothness tolerance is defined as the angular misalignment
    // at the points where the straight- and curved-line segments meet,
    // measured in radians. When computing the optimal path this quantity is
    // minimized, and the solver stops when reaching this tolerance.
    Real smoothnessTolerance = 0.1 / 180. * Pi; // Default = 0.1 degrees.
    // The solver's max allowed number of iterations used to compute the
    // optimal path.
    int solverMaxIterations = 50; // TODO a bit high?
    // The solver's max allowed stepsize measured in radians. This is converted
    // to a max allowed linear stepsize using the local radius of curvature
    // evaluated at each obstacle.
    Real solverMaxStepSize = 10. / 180. * Pi;
};

} // namespace

//==============================================================================
//                Contact Geometry Related Helper Functions
//==============================================================================
/* This section contains some helper functions for handling FrenetFrames and
flipping some signs in ContactGeometry. */
namespace
{

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

// Compute the frenet frame at a knot point of a geodesic.
FrenetFrame calcFrenetFrameFromGeodesicState(
    const ContactGeometry& geometry,
    const ContactGeometry::GeodesicKnotPoint& q)
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

// Compute the surface curvature such that tangentDot = curvature * normal
// (This flips the sign).
Real calcSurfaceCurvature(
    const ContactGeometry& geometry,
    const FrenetFrame& X_S,
    CoordinateAxis direction)
{
    return -geometry.calcSurfaceCurvatureInDirection(
        X_S.p(),
        X_S.R().getAxisUnitVec(direction));
}

// Evaluate whether a point lies above the given surface.
bool isPointAboveSurface(const ContactGeometry& geometry, const Vec3& point_S)
{
    // Try block was added here because some surfaces throw an exception when
    // calling calcSurfaceValue outside of a given range.
    try {
        // TODO surface value is negative above surface.
        return geometry.calcSurfaceValue(point_S) < 0.;
    } catch (...) {
        // TODO this is a bit too risky.
        return true;
    }
}

} // namespace

//==============================================================================
//                      Class CableSubsystem::Impl
//==============================================================================
// Implementation class of the CableSubsystem.
class CableSubsystem::Impl : public Subsystem::Guts {
public:
    //--------------------------------------------------------------------------
    //  Interface Translation Functions
    //--------------------------------------------------------------------------
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

    const CableSpan::Impl& getCableImpl(CableSpanIndex index) const;

    const MultibodySystem& getMultibodySystem() const
    {
        return MultibodySystem::downcast(getSystem());
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    //--------------------------------------------------------------------------
    //  Subsystem::Guts Virtual Interface
    //--------------------------------------------------------------------------
    Impl* cloneImpl() const override
    {
        return new Impl(*this);
    }

    int realizeSubsystemTopologyImpl(State& state) const override;

    int calcDecorativeGeometryAndAppendImpl(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const override;

    //--------------------------------------------------------------------------
    //  Data
    //--------------------------------------------------------------------------

    // All cables in this subsystem.
    Array_<CableSpan, CableSpanIndex> cables;

    //--------------------------------------------------------------------------
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                          Struct Curve Segment
//==============================================================================
/* The CableSpan's path consists of straight line segments between the
obstacles, and curved segments over the obstacles. This struct represents the
curve segment related to an obstacle in the path. It manages the cached data
related to that obstacle such as the WrappingStatus, the geodesic in local
coordinates, and the contact points in Ground frame coordinates. */
class CurveSegment final {
public:
    //--------------------------------------------------------------------------
    // Constructors.
    //--------------------------------------------------------------------------

    CurveSegment() = default;

    CurveSegment(
        CableSubsystem* subsystem,
        CableSpanIndex cableIndex,
        ObstacleIndex obstacleIndex,
        MobilizedBodyIndex obstacleBody,
        Transform X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry,
        const Vec3& contactPointHint_S) :
        m_subsystem(subsystem),
        m_cableIndex(cableIndex), m_obstacleIndex(obstacleIndex),
        m_body(obstacleBody), m_X_BS(std::move(X_BS)),
        m_geometry(std::move(obstacleGeometry)),
        m_contactPointHint_S(contactPointHint_S)
    {}

    //--------------------------------------------------------------------------
    // Realizing and cache access.
    //--------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& state)
    {
        CurveSegmentData::Instance dataInst;
        // Initialize the contact point if a hint is available.
        if (isContactPointHintAvailable()) {
            const Vec3 initContactPoint =
                getContactGeometry().projectDownhillToNearestPoint(
                    getContactPointHint());
            dataInst.X_SP.updP()    = initContactPoint;
            dataInst.X_SQ.updP()    = initContactPoint;
            dataInst.wrappingStatus = ObstacleWrappingStatus::InitialGuess;
        } else {
            // Otherwise assume no contact.
            dataInst.trackingPointOnLine_S = Vec3{0.};
            dataInst.wrappingStatus        = ObstacleWrappingStatus::LiftedFromSurface;
        }
        // Use auto-update discrete variable to retain the previous path as a
        // warmstart.
        m_indexDataInst = updSubsystem().allocateAutoUpdateDiscreteVariable(
            state,
            Stage::Position,
            new Value<CurveSegmentData::Instance>(dataInst),
            Stage::Instance);

        m_indexDataPos = updSubsystem().allocateCacheEntry(
            state,
            Stage::Position,
            Stage::Infinity,
            new Value<CurveSegmentData::Pos>());
    }

    void invalidatePosEntry(const State& state) const
    {
        getSubsystem().markCacheValueNotRealized(state, m_indexDataPos);
    }

    const CurveSegmentData::Instance& getDataInst(const State& state) const
    {
        const CableSubsystem& subsystem = getSubsystem();
        if (!subsystem.isDiscreteVarUpdateValueRealized(
                state,
                m_indexDataInst)) {
            updDataInst(state) = getPrevDataInst(state);
            subsystem.markDiscreteVarUpdateValueRealized(
                state,
                m_indexDataInst);
        }
        return Value<CurveSegmentData::Instance>::downcast(
            subsystem.getDiscreteVarUpdateValue(state, m_indexDataInst));
    }

    const CurveSegmentData::Pos& getDataPos(const State& state) const
    {
        return Value<CurveSegmentData::Pos>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataPos));
    }

    //------------------------------------------------------------------------------
    // Model configuration accessors.
    //------------------------------------------------------------------------------

    // Set the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    void setContactPointHint(Vec3 contactPointHint_S)
    {
        m_contactPointHint_S = contactPointHint_S;
    }

    // Get the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    Vec3 getContactPointHint() const
    {
        return m_contactPointHint_S;
    }

    bool isContactPointHintAvailable() const
    {
        return !(
            isNaN(m_contactPointHint_S[0]) && isNaN(m_contactPointHint_S[1]) &&
            isNaN(m_contactPointHint_S[2]));
    }

    const ContactGeometry& getContactGeometry() const
    {
        return *m_geometry;
    }

    void setContactGeometry(std::shared_ptr<const ContactGeometry> geometry)
    {
        m_geometry = std::move(geometry);
    }

    const MobilizedBody& getMobilizedBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_body);
    }

    const MobilizedBodyIndex& getMobilizedBodyIndex() const
    {
        return m_body;
    }

    void setMobilizedBodyIndex(MobilizedBodyIndex body)
    {
        m_body = body;
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
        if (!m_subsystem) {
            SimTK_ERRCHK_ALWAYS(
                m_subsystem,
                "CurveSegment::getSubsystem()",
                "CableSpan not yet adopted by any CableSubsystem");
        }
        return *m_subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        if (!m_subsystem) {
            SimTK_ERRCHK_ALWAYS(
                m_subsystem,
                "CurveSegment::updSubsystem()",
                "CableSpan not yet adopted by any CableSubsystem");
        }
        return *m_subsystem;
    }

    void setSubsystem(CableSubsystem& subsystem)
    {
        m_subsystem = &subsystem;
    }

    const CableSpan::Impl& getCable() const
    {
        return getSubsystem().getImpl().getCableImpl(m_cableIndex);
    }
    const IntegratorTolerances& getIntegratorTolerances() const;

    //------------------------------------------------------------------------------
    // Utility functions.
    //------------------------------------------------------------------------------

    Transform calcSurfaceFrameInGround(const State& state) const
    {
        return getMobilizedBody().getBodyTransform(state).compose(m_X_BS);
    }

    Vec3 calcInitialContactPoint_G(const State& state) const
    {
        return calcSurfaceFrameInGround(state).shiftFrameStationToBase(
            getDataInst(state).X_SP.p());
    }

    Vec3 calcFinalContactPoint_G(const State& state) const
    {
        return calcSurfaceFrameInGround(state).shiftFrameStationToBase(
            getDataInst(state).X_SQ.p());
    }

    bool isInContactWithSurface(const State& state) const
    {
        const CurveSegmentData::Instance& dataInst = getDataInst(state);
        return dataInst.wrappingStatus ==
                   ObstacleWrappingStatus::InContactWithSurface ||
               dataInst.wrappingStatus == ObstacleWrappingStatus::InitialGuess;
    }

    //------------------------------------------------------------------------------
    //  Helper Functions: Finding next and previous CurveSegments
    //------------------------------------------------------------------------------
    // The functions below essentially filter any obstacles that are not in
    // contact with the cable.

    // Find the first obstacle before this segment that is in contact with the
    // cable. Returns an invalid index if there is none.
    ObstacleIndex findPrevObstacleInContactWithCable(const State& state) const;
    // Similarly find the first obstacle after this segment that is in contact.
    ObstacleIndex findNextObstacleInContactWithCable(const State& state) const;

    // Find the last path point before this segment, in ground frame
    // coordinates.
    Vec3 findPrevPathPoint_G(const State& state) const;
    // Find the first path point after this segment, in ground frame
    // coordinates.
    Vec3 findNextPathPoint_G(const State& state) const;

    // Find the path point before this segment, in surface frame coordinates.
    Vec3 findPrevPathPoint_S(const State& state) const
    {
        const Vec3 prevPoint_S =
            calcSurfaceFrameInGround(state).shiftBaseStationToFrame(
                findPrevPathPoint_G(state));
        SimTK_ERRCHK_ALWAYS(
            isPointAboveSurface(getContactGeometry(), prevPoint_S),
            "CurveSegment::Impl::assertSurfaceBounds",
            "Preceding point lies inside the surface");
        return prevPoint_S;
    }
    // Find the path point after this segment, in surface frame coordinates.
    Vec3 findNextPathPoint_S(const State& state) const
    {
        const Vec3 nextPoint_S =
            calcSurfaceFrameInGround(state).shiftBaseStationToFrame(
                findNextPathPoint_G(state));
        SimTK_ERRCHK_ALWAYS(
            isPointAboveSurface(getContactGeometry(), nextPoint_S),
            "CurveSegment::Impl::assertSurfaceBounds",
            "Next point lies inside the surface");
        return nextPoint_S;
    }

    //------------------------------------------------------------------------------
    //  Updating CurveSegmentData::Instance Helpers
    //------------------------------------------------------------------------------

    // Compute a new geodesic from provided initial conditions.
    void shootNewGeodesic(
        const Vec3& point_S,
        const Vec3& tangent_S,
        Real length,
        Real initIntegratorStepSize,
        CurveSegmentData::Instance& dataInst) const
    {
        // Shoot a new geodesic over the surface, and compute the geodesic state
        // at the first and last contact point of the surface.
        std::array<ContactGeometry::GeodesicKnotPoint, 2>
            geodesicBoundaryStates;

        // First determine if this geometry has an analytic form.
        const ContactGeometry& geometry = getContactGeometry();
        const bool useAnalyticGeodesic  = geometry.isAnalyticFormAvailable();

        if (useAnalyticGeodesic) {
            // For analytic surfaces we can compute the inital and final frenet
            // frames without computing any intermediate knot points.
            int numberOfKnotPoints = 2;
            int knotIx             = 0;
            geometry.shootGeodesicInDirectionAnalytically(
                point_S,
                tangent_S,
                length,
                numberOfKnotPoints,
                // Copy the state at the boundary frames.
                [&](const ContactGeometry::GeodesicKnotPoint& q)
                {
                    geodesicBoundaryStates.at(knotIx) = q;
                    ++knotIx;
                });
        } else {
            // For implicit surfaces we must use a numerical integrator to
            // compute the initial and final frames. Because this is relatively
            // costly, we will cache this geodesic.
            dataInst.geodesicIntegratorStates.clear();
            const IntegratorTolerances& tols = getIntegratorTolerances();
            geometry.shootGeodesicInDirectionImplicitly(
                point_S,
                tangent_S,
                length,
                initIntegratorStepSize,
                tols.geodesicIntegratorAccuracy,
                tols.constraintProjectionTolerance,
                tols.constraintProjectionMaxIterations,
                [&](const ContactGeometry::GeodesicKnotPoint& q)
                {
                    // Store the knot points in the cache.
                    dataInst.geodesicIntegratorStates.push_back(q);
                });

            // Copy the state at the boundary frames.
            geodesicBoundaryStates = {
                dataInst.geodesicIntegratorStates.front(),
                dataInst.geodesicIntegratorStates.back(),
            };
        }

        // Use the geodeis state at the boundary frames to compute the fields in
        // CurveSegmentData::Instance.

        // Compute the Frenet frames at the boundary frames.
        dataInst.X_SP = calcFrenetFrameFromGeodesicState(
            geometry,
            geodesicBoundaryStates.front());
        dataInst.X_SQ = calcFrenetFrameFromGeodesicState(
            geometry,
            geodesicBoundaryStates.back());

        // Compute the jacobi field scalars at the final frame.
        dataInst.jacobi_Q = {
            geodesicBoundaryStates.back().jacobiTrans,
            geodesicBoundaryStates.back().jacobiRot,
        };
        dataInst.jacobiDot_Q = {
            geodesicBoundaryStates.back().jacobiTransDot,
            geodesicBoundaryStates.back().jacobiRotDot,
        };

        // Compute the curvatures at the boundary frames.
        dataInst.curvatures_P = {
            calcSurfaceCurvature(geometry, dataInst.X_SP, TangentAxis),
            calcSurfaceCurvature(geometry, dataInst.X_SP, BinormalAxis),
        };
        dataInst.curvatures_Q = {
            calcSurfaceCurvature(geometry, dataInst.X_SQ, TangentAxis),
            calcSurfaceCurvature(geometry, dataInst.X_SQ, BinormalAxis),
        };

        // Compute the geodesic torsion at the boundary frames.
        dataInst.torsion_P = geometry.calcSurfaceTorsionInDirection(
            dataInst.X_SP.p(),
            getTangent(dataInst.X_SP));
        dataInst.torsion_Q = geometry.calcSurfaceTorsionInDirection(
            dataInst.X_SQ.p(),
            getTangent(dataInst.X_SQ));

        // Store the arc length.
        dataInst.arcLength = length;

        // Update the initial integrator step size:
        if (!useAnalyticGeodesic) { // does not apply to analytic surfaces.

            // We want to know if we should attempt a larger initial step size
            // next time, or a smaller one. If the initIntegratorStepSize was
            // rejected, it makes little sense to start with the exact same step
            // size next time, only to find it rejected again.

            const int numStepsTaken =
                static_cast<int>(dataInst.geodesicIntegratorStates.size()) - 1;

            // If the integrator took a single step or none: The final arc
            // length is shorter than the initial step size attempted. The step
            // was atleast not rejected, but we do not know if we can make it
            // larger. So we do not update it.
            if (numStepsTaken <= 1) {
                // Try the same step next time:
                dataInst.integratorInitialStepSize = initIntegratorStepSize;
            }

            // If the integrator took two or more steps, we might want to update
            // the next step size to try,
            if (numStepsTaken >= 2) {
                // The acual initial step size taken by the integrator cannot
                // have been larger than initIntegratorStepSize, but it might be
                // smaller if it was rejected. If rejected, we will reduce the
                // init step size for next time.
                const Real actualFirstStepSize =
                    dataInst.geodesicIntegratorStates.at(1).arcLength -
                    dataInst.geodesicIntegratorStates.front().arcLength;
                const bool initStepWasRejected =
                    actualFirstStepSize < initIntegratorStepSize;
                if (initStepWasRejected) {
                    // Reduce the init step size for next time, to avoid
                    // rejecting it again.
                    dataInst.integratorInitialStepSize = actualFirstStepSize;
                } else {
                    // If the initIntegratorStepSize was accepted we might want
                    // to try a larger step next time. The actualFirstStepSize
                    // cannot be larger than initIntegratorStepSize: we must
                    // take a look at the second step size.
                    const Real actualSecondStepSize =
                        dataInst.geodesicIntegratorStates.at(2).arcLength -
                        dataInst.geodesicIntegratorStates.at(1).arcLength;
                    // We already established that the initIntegratorStepSize
                    // was accepted, now we are only looking to increase the
                    // step size.
                    dataInst.integratorInitialStepSize =
                        std::max(initIntegratorStepSize, actualSecondStepSize);
                }
            }
        }

        dataInst.wrappingStatus = ObstacleWrappingStatus::InContactWithSurface;
    }

    // Lift curve from surface, and start tracking the given point.
    void liftCurveFromSurface(
        const Vec3& trackingPoint_S,
        CurveSegmentData::Instance& dataInst) const
    {
        dataInst.wrappingStatus = ObstacleWrappingStatus::LiftedFromSurface;
        dataInst.arcLength      = NaN;
        dataInst.trackingPointOnLine_S = trackingPoint_S;
    }

    // Helper function for initializing the path.
    void calcInitCurve(const State& state) const
    {
        // Check if a valid contact point hint was set.
        if (!isContactPointHintAvailable()) {
            // No initial contact point was set: Initialize as lifted from the
            // surface, and initialize tracking point at origin (for lack of
            // alternative).
            liftCurveFromSurface(Vec3(0.), updDataInst(state));
            return;
        }
        // If we do have a contact point hint, initialize the path by shooting a
        // zero-length geodesic at that location.

        // Shoot a zero length geodesic at the contact point with tangent
        // directed from prevPoint to nextPoint:
        const Vec3 prevPoint_S = findPrevPathPoint_S(state);
        const Vec3 nextPoint_S = findNextPathPoint_S(state);
        shootNewGeodesic(
            getContactPointHint(),     // = location
            nextPoint_S - prevPoint_S, // = tangent
            0.,                        // = arc length
            getDataInst(state).integratorInitialStepSize,
            updDataInst(state));
    }

    // Helper function for detecting touchdown on the obstacle's surface, and
    // the location of the touchdown point if so.
    void calcCurveTouchdownIfNeeded(const State& state) const
    {
        // Only attempt touchdown when lifted.
        const CurveSegmentData::Instance& dataInst = getDataInst(state);
        if (dataInst.wrappingStatus !=
            ObstacleWrappingStatus::LiftedFromSurface) {
            return;
        }

        // Given the straight line between the point before and after this
        // segment, touchdown is detected by computing the point on that line
        // that is nearest to the surface.
        const Vec3 prevPoint_S = findPrevPathPoint_S(state);
        const Vec3 nextPoint_S = findNextPathPoint_S(state);

        // Use the cached tracking point as the initial guess.
        Vec3 pointOnLineNearSurface_S = dataInst.trackingPointOnLine_S;

        // Compute the point on the line nearest the surface.
        const IntegratorTolerances& tols = getIntegratorTolerances();
        bool touchdownDetected           = false;
        try {
            // This routine might fail to converge if the line is too far from
            // the surface. This should not be a problem, because in that case
            // we definitely have not touched down.
            touchdownDetected =
                getContactGeometry().calcNearestPointOnLineImplicitly(
                    prevPoint_S,
                    nextPoint_S,
                    tols.constraintProjectionMaxIterations,
                    tols.constraintProjectionTolerance,
                    pointOnLineNearSurface_S);
        } catch (...) {
        }

        // In case of touchdown, shoot a zero-length geodesic at the touchdown
        // point.
        if (touchdownDetected) {
            shootNewGeodesic(
                pointOnLineNearSurface_S,
                nextPoint_S - prevPoint_S,
                0.,
                dataInst.integratorInitialStepSize,
                updDataInst(state));
            return;
        }

        // Not touchingdown indicates liftoff:
        liftCurveFromSurface(pointOnLineNearSurface_S, updDataInst(state));
    }

    // Helper function for detecting liftoff from the obstacle's surface.
    void calcCurveLiftoffIfNeeded(const State& state) const
    {
        // Only attempt liftoff when currently wrapping the surface.
        const CurveSegmentData::Instance& dataInst = getDataInst(state);
        if (dataInst.wrappingStatus !=
            ObstacleWrappingStatus::InContactWithSurface) {
            return;
        }

        // The curve length must have shrunk completely before lifting off.
        if (dataInst.arcLength > 0.) {
            return;
        }

        // For a zero-length curve, trigger liftoff when the prev and next
        // points lie above the surface plane.
        const Vec3 prevPoint_S = findPrevPathPoint_S(state);
        const Vec3 nextPoint_S = findNextPathPoint_S(state);
        if (dot(prevPoint_S - dataInst.X_SP.p(),
                dataInst.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0. ||
            dot(nextPoint_S - dataInst.X_SP.p(),
                dataInst.X_SP.R().getAxisUnitVec(NormalAxis)) <= 0.) {
            // No liftoff.
            return;
        }

        // Liftoff detected, initialize the tracking point from the last contact
        // point.
        liftCurveFromSurface(dataInst.X_SP.p(), updDataInst(state));
    }

    //------------------------------------------------------------------------------
    //  Computing cached data
    //------------------------------------------------------------------------------

    // Apply the correction to the initial condition of the geodesic, and
    // shoot a new geodesic, updating the cache variable.
    void applyGeodesicCorrection(
        const State& state,
        const NaturalGeodesicCorrection& c) const
    {
        // Get the previous geodesic.
        const CurveSegmentData::Instance& dataInst = getDataInst(state);

        // Frenet frame at initial contact point.
        const FrenetFrame& X_SP = dataInst.X_SP;
        const UnitVec3& t       = getTangent(X_SP);
        const UnitVec3& n       = getNormal(X_SP);
        const UnitVec3& b       = getBinormal(X_SP);

        const Real tau   = dataInst.torsion_P;
        const Real kappa = dataInst.curvatures_P[0];

        // Get corrected initial conditions.
        const Vec3 dx               = t * c[0] + b * c[1];
        const Vec3 correctedPoint_S = X_SP.p() + dx;

        const Vec3 w = (kappa * c[0] - tau * c[1]) * b - c[2] * n;
        const Vec3 correctedTangent_S = t + cross(w, t);

        // Take the length correction, and add to the current length.
        const Real dl = c[3]; // Length increment is the last element.
        const Real correctedLength = std::max(
            dataInst.arcLength + dl,
            0.); // Clamp length to be nonnegative.

        // Shoot the new geodesic.
        shootNewGeodesic(
            correctedPoint_S,
            correctedTangent_S,
            correctedLength,
            dataInst.integratorInitialStepSize,
            updDataInst(state));

        getSubsystem().markDiscreteVarUpdateValueRealized(
            state,
            m_indexDataInst);
        invalidatePosEntry(state);
    }

    const CurveSegmentData::Instance& calcDataInst(const State& state) const
    {
        // Update the CurveSegmentData::Instance cache depending on the current
        // wrapping status.
        switch (getDataInst(state).wrappingStatus) {
        case ObstacleWrappingStatus::Disabled:
            break;
        case ObstacleWrappingStatus::InitialGuess:
            // Path needs to be initialized.
            calcInitCurve(state);
            break;
        case ObstacleWrappingStatus::InContactWithSurface:
            // If in contact with obstacle: detect lifting off.
            calcCurveLiftoffIfNeeded(state);
            break;
        case ObstacleWrappingStatus::LiftedFromSurface:
            // If not in contact with obstacle: detect touching down.
            calcCurveTouchdownIfNeeded(state);
            break;
        }

        getSubsystem().markDiscreteVarUpdateValueRealized(
            state,
            m_indexDataInst);
        return getDataInst(state);
    }

    const CurveSegmentData::Pos& calcDataPos(const State& state) const
    {
        // Make sure the Instance level data is up to date.
        const CurveSegmentData::Instance& dataInst = calcDataInst(state);

        // Update the Stage::Position level cache.
        CurveSegmentData::Pos& dataPos = updDataPos(state);

        // Check if the obstacle is enabled.
        if (dataInst.wrappingStatus != ObstacleWrappingStatus::Disabled) {
            // Store tramsform from local surface frame to ground.
            dataPos.X_GS = calcSurfaceFrameInGround(state);
            // Store the geodesic's frenet frames in ground frame.
            dataPos.X_GP = dataPos.X_GS.compose(dataInst.X_SP);
            dataPos.X_GQ = dataPos.X_GS.compose(dataInst.X_SQ);
        }

        getSubsystem().markCacheValueRealized(state, m_indexDataPos);

        return dataPos;
    }

    //------------------------------------------------------------------------------
    //  Generating Geodesic Points For Visualization
    //------------------------------------------------------------------------------

    void calcDecorativePathPoints(
        const State& state,
        const std::function<void(Vec3 point_G)>& sink) const
    {
        if (!isInContactWithSurface(state)) {
            return;
        }

        const CurveSegmentData::Instance& dataInst = getDataInst(state);
        const Transform& X_GS                      = getDataPos(state).X_GS;

        // Lambda for computing the point in ground from a GeodesicKnotPoint,
        // and writing it to the sink.
        auto calcPathPointAndWriteToSink =
            [&](const ContactGeometry::GeodesicKnotPoint& q)
        {
            const Vec3 point_G = X_GS.shiftFrameStationToBase(q.point);
            sink(point_G);
        };

        if (getContactGeometry().isAnalyticFormAvailable()) {
            // TODO Depending on what is actually desired by the end-user, this
            // might need revision.

            // For analytic surfaces there is no numerical integrator to give a
            // natural interspacing of knots. So instead I just assume you
            // would like some OK granularity. This can be done by interspacing
            // the knots by a certain "angle" (e.g. 5 degrees is pretty fine),
            // and using the radius of curvature to convert the angle to a
            // spatial interspacing of the samples.
            constexpr Real c_AngularSpacingInDegrees = 5.;

            // To convert the angular spacing to a spatial spacing we need the
            // curvature. The curvature is generally not constant along a
            // geodesic, but there are not a lot of analytic surfaces at the
            // moment (sphere and cylinder?). For these we can just look at the
            // curvature at the boundary frames to know the curvature along the
            // entire geodesic.
            const Real maxCurvature_P = max(abs(dataInst.curvatures_P));
            const Real maxCurvature_Q = max(abs(dataInst.curvatures_P));
            const Real minRadiusOfCurvature =
                1. / std::max(maxCurvature_P, maxCurvature_Q);

            // Estimate the spatial interspacing from the angular spacing and
            // the radius of curvature.
            const Real lengthInterspacing =
                minRadiusOfCurvature * c_AngularSpacingInDegrees / 180. * Pi;

            // Now shoot the geodesic analytically with the desired granularity.
            const int numberOfKnotPoints =
                static_cast<int>(dataInst.arcLength / lengthInterspacing);
            getContactGeometry().shootGeodesicInDirectionAnalytically(
                dataInst.X_SP.p(),
                getTangent(dataInst.X_SP),
                dataInst.arcLength,
                std::max(numberOfKnotPoints, 2),
                calcPathPointAndWriteToSink);
            return;
        }

        // If the surface did not have an analytic form, we simply forward the
        // knots from the GeodesicIntegrator.
        for (const ContactGeometry::GeodesicKnotPoint& q :
             dataInst.geodesicIntegratorStates) {
            calcPathPointAndWriteToSink(q);
        }
    }

private:
    //------------------------------------------------------------------------------
    // Private Cache Access.
    //------------------------------------------------------------------------------

    CurveSegmentData::Instance& updDataInst(const State& state) const
    {
        return Value<CurveSegmentData::Instance>::updDowncast(
            getSubsystem().updDiscreteVarUpdateValue(state, m_indexDataInst));
    }

    const CurveSegmentData::Instance& getPrevDataInst(const State& state) const
    {
        return Value<CurveSegmentData::Instance>::downcast(
            getSubsystem().getDiscreteVariable(state, m_indexDataInst));
    }

    CurveSegmentData::Instance& updPrevDataInst(State& state) const
    {
        return Value<CurveSegmentData::Instance>::updDowncast(
            getSubsystem().updDiscreteVariable(state, m_indexDataInst));
    }

    CurveSegmentData::Pos& updDataPos(const State& state) const
    {
        return Value<CurveSegmentData::Pos>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataPos));
    }

    //------------------------------------------------------------------------------
    // Data
    //------------------------------------------------------------------------------

    // Subsystem info.
    CableSubsystem* m_subsystem = nullptr;
    // The index of this CurveSegment's cable in the CableSubsystem.
    CableSpanIndex m_cableIndex = CableSpanIndex::Invalid();
    // The index of this CurveSegment's obstacle in the cable.
    ObstacleIndex m_obstacleIndex = ObstacleIndex::Invalid();

    // MobilizedBody that surface is attached to.
    MobilizedBodyIndex m_body;
    // Surface to body transform.
    Transform m_X_BS;

    // Obstacle surface.
    std::shared_ptr<const ContactGeometry> m_geometry;

    // Topology cache.
    CacheEntryIndex m_indexDataPos        = CacheEntryIndex::Invalid();
    DiscreteVariableIndex m_indexDataInst = DiscreteVariableIndex::Invalid();

    // Initial contact point hint used to setup the initial path.
    Vec3 m_contactPointHint_S{NaN, NaN, NaN};

    //------------------------------------------------------------------------------
    // Helper class for unit tests.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                         Class Cablespan::Impl
//==============================================================================
// TODO describe (the impl side of CableSpan)
class CableSpan::Impl {
public:
    //--------------------------------------------------------------------------
    // Constructors.
    //--------------------------------------------------------------------------

    Impl() = default;

    // Construct a CableSpan::Impl with valid data fields, but without being
    // registered at any CableSubsystem.
    Impl(
        MobilizedBodyIndex originBody,
        Vec3 originPoint_B,
        MobilizedBodyIndex terminationBody,
        Vec3 terminationPoint_B) :
        m_originBody(originBody),
        m_originPoint_B(originPoint_B), m_terminationBody(terminationBody),
        m_terminationPoint_B(terminationPoint_B)
    {}

    //--------------------------------------------------------------------------
    // Realizing and accessing cache.
    //--------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& state, CacheEntryIndex indexOfDataInstEntry)
    {
        m_indexDataInst = indexOfDataInstEntry;

        m_indexDataPos = updSubsystem().allocateCacheEntry(
            state,
            Stage::Position,
            Stage::Infinity,
            new Value<CableSpanData::Pos>());

        m_indexDataVel = updSubsystem().allocateCacheEntry(
            state,
            Stage::Velocity,
            Stage::Infinity,
            new Value<CableSpanData::Vel>());

        for (CurveSegment& segment : m_curveSegments) {
            segment.realizeTopology(state);
        }
    }

    void invalidateTopology()
    {
        if (m_subsystem) {
            getSubsystem().invalidateSubsystemTopologyCache();
        }
    }

    void realizePosition(const State& state) const
    {
        if (getSubsystem().isCacheValueRealized(state, m_indexDataPos)) {
            return;
        }
        calcDataPos(state);
        getSubsystem().markCacheValueRealized(state, m_indexDataPos);
    }

    void realizeVelocity(const State& state) const
    {
        realizePosition(state);
        if (getSubsystem().isCacheValueRealized(state, m_indexDataVel)) {
            return;
        }
        calcDataVel(state, updDataVel(state));
        getSubsystem().markCacheValueRealized(state, m_indexDataVel);
    }

    const CableSpanData::Pos& getDataPos(const State& state) const
    {
        realizePosition(state);
        return Value<CableSpanData::Pos>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataPos));
    }

    const CableSpanData::Vel& getDataVel(const State& state) const
    {
        realizeVelocity(state);
        return Value<CableSpanData::Vel>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataVel));
    }

    //--------------------------------------------------------------------------
    // Accessors
    //--------------------------------------------------------------------------

    const CurveSegment& getObstacleCurveSegment(ObstacleIndex ix) const
    {
        return m_curveSegments[ix];
    }

    CurveSegment& updObstacleCurveSegment(ObstacleIndex ix)
    {
        return m_curveSegments[ix];
    }

    // Get the origin point in body fixed coordinates.
    const Vec3& getOriginPoint_B() const
    {
        return m_originPoint_B;
    }

    // Get the termination point in body fixed coordinates.
    const Vec3& getTerminationPoint_B() const
    {
        return m_terminationPoint_B;
    }

    CableSpanParameters& updParameters()
    {
        return m_parameters;
    }

    const CableSpanParameters& getParameters() const
    {
        return m_parameters;
    }

    //--------------------------------------------------------------------------
    // Helper functions: End-points frame transformations.
    //--------------------------------------------------------------------------

    // Get body to which this cable's origin point is rigidly attached to.
    const Mobod& getOriginBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_originBody);
    }

    // Get body to which this cable's termination point is rigidly attached to.
    const Mobod& getTerminationBody() const
    {
        return getSubsystem()
            .getMultibodySystem()
            .getMatterSubsystem()
            .getMobilizedBody(m_terminationBody);
    }

    Vec3 calcOriginPointInGround(const State& state) const
    {
        return getOriginBody().getBodyTransform(state).shiftFrameStationToBase(
            m_originPoint_B);
    }

    Vec3 calcTerminationPointInGround(const State& state) const
    {
        return getTerminationBody()
            .getBodyTransform(state)
            .shiftFrameStationToBase(m_terminationPoint_B);
    }

    //--------------------------------------------------------------------------
    // Helpers functions: Accessing obstacles in contact with cable.
    //--------------------------------------------------------------------------

    // Iterate over the cable's curve segments, and call provided function for
    // those that are in contact with the obstacle's surfacce.
    void forEachActiveCurveSegment(
        const State& s,
        const std::function<void(const CurveSegment& curve)>& callMe) const
    {
        for (ObstacleIndex ix(0); ix < getNumObstacles(); ++ix) {
            const CurveSegment& curve = getObstacleCurveSegment(ix);
            if (!curve.isInContactWithSurface(s)) {
                continue;
            }
            callMe(curve);
        }
    }

    // Count the number of CurveSegments that are in contact with the obstacle's
    // surface.
    int countActive(const State& s) const
    {
        int counter = 0;
        forEachActiveCurveSegment(s, [&](const CurveSegment&) { ++counter; });
        return counter;
    }

    //--------------------------------------------------------------------------
    // CableSpan interface translation.
    //--------------------------------------------------------------------------

    ObstacleIndex addObstacle(
        MobilizedBodyIndex obstacleBody,
        const Transform& X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry,
        const Vec3& contactPointHint_S)
    {
        invalidateTopology();

        ObstacleIndex obstacleIx(m_curveSegments.size());

        m_curveSegments.push_back(CurveSegment(
            m_subsystem,
            getIndex(),
            obstacleIx,
            obstacleBody,
            X_BS,
            std::move(obstacleGeometry),
            contactPointHint_S));

        return obstacleIx;
    }

    int getNumObstacles() const
    {
        return m_curveSegments.size();
    }

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    Real calcCablePower(const State& state, Real tension) const;

    void calcDecorativePathPoints(
        const State& state,
        const std::function<void(Vec3 point_G)>& sink) const
    {
        // Write the initial path point.
        const CableSpanData::Pos& dataPos = getDataPos(state);
        sink(dataPos.originPoint_G);

        // Write path points along each of the curves.
        for (const CurveSegment& curve : m_curveSegments) {
            curve.calcDecorativePathPoints(state, sink);
        }

        // Write the path's termination point.
        sink(dataPos.terminationPoint_G);
    }

    void calcDecorativeGeometryAndAppend(
        const State& state,
        Array_<DecorativeGeometry>& decorations) const
    {
        // Draw lines between all path points.
        {
            const Vec3 c_Color            = Purple;
            constexpr int c_LineThickness = 3;

            bool isFirstPoint = true;
            Vec3 prevPoint_G{NaN};
            // Lambda that draws a line between two path points prevPoint_G &
            // nextPoint_G and pushes it to the decorations buffer.
            auto pushDecorativeLineBetweenPathPoints = [&](Vec3 nextPoint_G)
            {
                if (!isFirstPoint) {
                    decorations.push_back(
                        DecorativeLine(prevPoint_G, nextPoint_G)
                            .setColor(c_Color)
                            .setLineThickness(c_LineThickness));
                }
                prevPoint_G  = nextPoint_G;
                isFirstPoint = false;
            };
            // Call the lambda for all path points.
            calcDecorativePathPoints(
                state,
                pushDecorativeLineBetweenPathPoints);
        }

        // For any obstacle that is not in contact with the cable, draw the
        // tracking point used for touchdown.
        {
            const Vec3 c_Color            = Red;
            constexpr int c_LineThickness = 3;

            // Find those obstacles that are not in contact.
            for (ObstacleIndex ix(0); ix < getNumObstacles(); ++ix) {
                const CurveSegment& curve = getObstacleCurveSegment(ix);
                if (curve.isInContactWithSurface(state)) {
                    continue;
                }

                // Draw the tracking point (in ground frame).
                const Transform X_GS = curve.getDataPos(state).X_GS;
                const Vec3 trackingPoint_S =
                    curve.getDataInst(state).trackingPointOnLine_S;
                const Vec3 trackingPoint_G =
                    X_GS.shiftFrameStationToBase(trackingPoint_S);
                decorations.push_back(
                    DecorativePoint(trackingPoint_G).setColor(c_Color));
                // Draw line to tracking point from either the surface, or the
                // obstacle origin. We default to drawing from the origin
                // because surface projection becomes unstable at larger
                // distances, which causes the projection to take too much
                // time. Which was the reason to not use it in the first place.
                // But it is more illustrative.
                constexpr bool c_DrawLineFromObstacleOrigin = true;
                if (c_DrawLineFromObstacleOrigin) {
                    decorations.push_back(
                        DecorativeLine(X_GS.p(), trackingPoint_G)
                            .setColor(c_Color)
                            .setLineThickness(c_LineThickness));
                } else {
                    try {
                        const Vec3 surfacePoint_S =
                            curve.getContactGeometry()
                                .projectDownhillToNearestPoint(trackingPoint_S);
                        const Vec3 surfacePoint_G =
                            X_GS.shiftFrameStationToBase(surfacePoint_S);
                        decorations.push_back(
                            DecorativeLine(trackingPoint_G, surfacePoint_G)
                                .setColor(Orange)
                                .setLineThickness(c_LineThickness));
                    } catch (...) {
                    }
                }
            }
        }
    }

private:
    //--------------------------------------------------------------------------
    // Private accessors.
    //--------------------------------------------------------------------------

    // Get the subsystem this CableSpan is part of.
    const CableSubsystem& getSubsystem() const
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CableSpan::Impl::getSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");
        return *m_subsystem;
    }

    void setSubsystem(CableSubsystem& subsystem, CableSpanIndex index)
    {
        m_subsystem = &subsystem;
        m_index     = index;
        for (CurveSegment& curve : m_curveSegments) {
            curve.setSubsystem(subsystem);
        }
    }

    CableSubsystem& updSubsystem()
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CableSpan::Impl::updSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");
        return *m_subsystem;
    }

    // Get the index of this CableSpan in the CableSubsystem.
    CableSpanIndex getIndex() const
    {
        SimTK_ERRCHK1_ALWAYS(
            m_index.isValid(),
            "CableSpan::Impl::getIndex()",
            "Index %d is not valid.",
            m_index);
        return m_index;
    }

    // Mutable CableSpanData::Instance cache access.
    CableSpanData::Instance& updDataInst(const State& state) const
    {
        return Value<CableSpanData::Instance>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataInst));
    }

    // Mutable CableSpanData::Pos cache access.
    CableSpanData::Pos& updDataPos(const State& state) const
    {
        return Value<CableSpanData::Pos>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataPos));
    }

    // Mutable CableSpanData::Vel cache access.
    CableSpanData::Vel& updDataVel(const State& state) const
    {
        return Value<CableSpanData::Vel>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataVel));
    }

    //------------------------------------------------------------------------------
    //  Computing cached data
    //------------------------------------------------------------------------------

    const MatrixWorkspace& calcDataInst(const State& state) const;

    const CableSpanData::Pos& calcDataPos(const State& s) const
    {
        CableSpanData::Pos& dataPos = updDataPos(s);

        dataPos.originPoint_G      = calcOriginPointInGround(s);
        dataPos.terminationPoint_G = calcTerminationPointInGround(s);

        // Helper for computing the total cable length.
        auto calcTotalCableLength =
            [&](const std::vector<LineSegment>& lines) -> Real
        {
            Real totalCableLength = 0.;
            // Add length of each line segment.
            for (const LineSegment& line : lines) {
                totalCableLength += line.length;
            }
            // Add length of each curve segment.
            forEachActiveCurveSegment(
                s,
                [&](const CurveSegment& curve)
                {
                    // Update the total cable length.
                    totalCableLength += curve.getDataInst(s).arcLength;
                });
            return totalCableLength;
        };

        // Computing the CableSpan's path is done iteratively by computing
        // corrections for the CurveSegments until the pathErrorVector is small
        // enough.
        dataPos.loopIter = 0;
        while (true) {
            // Compute all data required for updating the CableSpanData::Pos
            // cache.
            const MatrixWorkspace& data = calcDataInst(s);

            // Stop iterating if:
            // - No obstacles in contact with path: Path is straight line.
            // - Converged: No significant correction.
            // - Not converged: Max iterations has been reached.
            if (data.nObstaclesInContact == 0 ||
                data.pathCorrection.normInf() < Eps ||
                dataPos.loopIter >= getParameters().solverMaxIterations) {
                // Update cache entry and stop solver.
                dataPos.smoothness  = data.maxPathError;
                dataPos.cableLength = calcTotalCableLength(data.lineSegments);
                dataPos.originTangent_G = data.lineSegments.front().direction;
                dataPos.terminationTangent_G =
                    data.lineSegments.back().direction;
                break;
            }

            // Apply the corrections to each CurveSegment to reduce the path
            // error.
            int activeCurveIx = 0; // Index of curve in all that are in contact.
            forEachActiveCurveSegment(
                s,
                [&](const CurveSegment& curve)
                {
                    curve.applyGeodesicCorrection(
                        s,
                        data.getCurveCorrection(activeCurveIx));
                    ++activeCurveIx;
                });

            // The applied corrections have changed the path: invalidate each
            // segment's cache.
            for (const CurveSegment& curve : m_curveSegments) {
                // Also invalidate non-active segments: They might touchdown
                // again.
                curve.invalidatePosEntry(s);
            }

            ++dataPos.loopIter;
        }

        return dataPos;
    }

    void calcDataVel(const State& s, CableSpanData::Vel& dataVel) const
    {
        const CableSpanData::Pos& dataPos = getDataPos(s);

        auto CalcPointVelocityInGround = [&](const MobilizedBody& mobod,
                                             const Vec3& point_G) -> Vec3
        {
            // Not using MobilizedBody::findStationVelocityInGround because the
            // point is already in ground frame.

            // Get body kinematics in ground frame.
            const Vec3& x_BG = mobod.getBodyOriginLocation(s);
            const Vec3& w_BG = mobod.getBodyAngularVelocity(s);
            const Vec3& v_BG = mobod.getBodyOriginVelocity(s);

            // Compute surface point velocity in ground frame.
            return v_BG + w_BG % (point_G - x_BG);
        };

        Real& lengthDot = (dataVel.lengthDot = 0.);

        Vec3 v_GQ =
            getOriginBody().findStationVelocityInGround(s, m_originPoint_B);
        Vec3 x_GQ = m_originPoint_B;

        forEachActiveCurveSegment(
            s,
            [&](const CurveSegment& curve)
            {
                const MobilizedBody& mobod = curve.getMobilizedBody();
                const CurveSegmentData::Pos& curveDataPos = curve.getDataPos(s);

                const UnitVec3& e_G = getTangent(curveDataPos.X_GP);

                const Vec3 v_GP =
                    CalcPointVelocityInGround(mobod, curveDataPos.X_GP.p());

                lengthDot += dot(e_G, v_GP - v_GQ);

                x_GQ = curveDataPos.X_GQ.p();
                v_GQ = CalcPointVelocityInGround(mobod, x_GQ);
            });

        const Vec3 v_GP = getTerminationBody().findStationVelocityInGround(
            s,
            m_terminationPoint_B);

        const UnitVec3 e_G(dataPos.terminationPoint_G - x_GQ);

        lengthDot += dot(e_G, v_GP - v_GQ);
    }

    //------------------------------------------------------------------------------
    // Data
    //------------------------------------------------------------------------------

    // Reference back to the subsystem.
    CableSubsystem* m_subsystem = nullptr;
    CableSpanIndex m_index      = CableSpanIndex::Invalid();

    // TOPOLOGY CACHE (set during realizeTopology())
    CacheEntryIndex m_indexDataInst = CacheEntryIndex::Invalid();
    CacheEntryIndex m_indexDataPos  = CacheEntryIndex::Invalid();
    CacheEntryIndex m_indexDataVel  = CacheEntryIndex::Invalid();

    // Path origin point (body and point in body coordinates).
    MobilizedBodyIndex m_originBody = MobilizedBodyIndex::Invalid();
    Vec3 m_originPoint_B{NaN};

    // Path termination point (body and point in body coordinates).
    MobilizedBodyIndex m_terminationBody = MobilizedBodyIndex::Invalid();
    Vec3 m_terminationPoint_B{NaN};

    // Path CurveSegments over the obstacles.
    Array_<CurveSegment, ObstacleIndex> m_curveSegments{};

    // All configuration parameters.
    CableSpanParameters m_parameters{};

    friend CableSpan;
    friend CableSubsystem;
    friend CableSubsystemTestHelper;
};

//==============================================================================
//          CurveSegment Utility Functions Using CableSpan Handle
//==============================================================================
// These helper functions are defined here as they depend on the declarations
// in CableSpan::Impl.

const IntegratorTolerances& CurveSegment::getIntegratorTolerances() const
{
    return getCable().getParameters();
}

ObstacleIndex CurveSegment::findPrevObstacleInContactWithCable(
    const State& state) const
{
    const CableSpan::Impl& cable = getCable();
    // Find the first obstacle that makes contact before this CurveSegment.
    for (ObstacleIndex ix(m_obstacleIndex); ix > 0;) {
        --ix;
        if (cable.getObstacleCurveSegment(ix).isInContactWithSurface(state)) {
            return ix;
        }
    }
    // No obstacles in contact before this segment.
    return ObstacleIndex::Invalid();
}

ObstacleIndex CurveSegment::findNextObstacleInContactWithCable(
    const State& state) const
{
    const CableSpan::Impl& cable = getCable();
    // Find the first obstacle that makes contact after this CurveSegment.
    for (ObstacleIndex ix(m_obstacleIndex + 1); ix < cable.getNumObstacles();
         ++ix) {
        if (cable.getObstacleCurveSegment(ix).isInContactWithSurface(state)) {
            return ix;
        }
    }
    // No obstacles in contact after this segment.
    return ObstacleIndex::Invalid();
}

Vec3 CurveSegment::findPrevPathPoint_G(const State& state) const
{
    // Check if there is a curve segment preceding given obstacle.
    const ObstacleIndex prevObstacle =
        findPrevObstacleInContactWithCable(state);
    if (prevObstacle.isValid()) {
        // The previous point is the final contact point of the previous curve.
        return getCable()
            .getObstacleCurveSegment(prevObstacle)
            .calcFinalContactPoint_G(state);
    }
    // There are no curve segments before given obstacle: the previous point
    // is the path's origin point.
    const CableSpan::Impl& cable = getCable();
    return getCable()
        .getOriginBody()
        .getBodyTransform(state)
        .shiftFrameStationToBase(cable.getOriginPoint_B());
}

Vec3 CurveSegment::findNextPathPoint_G(const State& state) const
{
    // Check if there is a curve segment after given obstacle.
    const ObstacleIndex nextObstacle =
        findNextObstacleInContactWithCable(state);
    if (nextObstacle.isValid()) {
        // The next point is the initial contact point of the next curve.
        return getCable()
            .getObstacleCurveSegment(nextObstacle)
            .calcInitialContactPoint_G(state);
    };
    // There are no curve segments following given obstacle: the next point
    // is the path's termination point.
    const CableSpan::Impl& cable = getCable();
    return cable.getTerminationBody()
        .getBodyTransform(state)
        .shiftFrameStationToBase(cable.getTerminationPoint_B());
}

//==============================================================================
//                      CABLE FORCE COMPUTATIONS
//==============================================================================
// This section contains all force related computations.
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
    const CableSpanData::Pos& dataPos = cable.getDataPos(s);

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
    const CableSpanData::Pos& dataPos = cable.getDataPos(s);

    const Vec3& arm_G = dataPos.terminationPoint_G -
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
    for (const CurveSegment& curve : m_curveSegments) {
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

    for (const CurveSegment& curve : m_curveSegments) {
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
//                             Path Error Vector
//==============================================================================
// This section contains helper functions for computing the path error vector.
namespace
{

// Compute the straight line segments in between the obstacles.
void calcLineSegments(
    const CableSpan::Impl& cable,
    const State& s,
    const Vec3& pathOriginPoint,
    const Vec3& pathTerminationPoint,
    std::vector<LineSegment>& lines)
{
    lines.clear();
    Vec3 prevPathPoint(pathOriginPoint);

    cable.forEachActiveCurveSegment(
        s,
        [&](const CurveSegment& curve)
        {
            const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);

            // Compute the line segment from the previous point to the
            // first curve contact point.
            lines.emplace_back(prevPathPoint, dataPos.X_GP.p());

            // The next line segment will start at the curve's final contact
            // point.
            prevPathPoint = dataPos.X_GQ.p();
        });

    // Compute the last line segment.
    lines.emplace_back(prevPathPoint, pathTerminationPoint);
}

// Helper for computing a single element of the path error vector.
// The algorithm for finding the correct path will attempt to drive this path
// error to zero. Since we would like the curve's tangents to be aligned with
// the straight line segments, the path error is defined as the alignment of the
// line segments with the normal and binormal axes.
// In that case: zero path error -> no alignment of line segments with normal
// and binormal -> line segments must be aligned with the tangents.
Real calcPathErrorElement(
    const LineSegment& e,
    const FrenetFrame& X,
    CoordinateAxis axis)
{
    return dot(e.direction, X.R().getAxisUnitVec(axis));
}

// Compute the path error for each curve segment at the boundary frames, and
// stack them in a vector. At each boundary frame (X_GP, X_GQ) the path error
// along the given axes is computed, in that order.
template <size_t N>
void calcPathErrorVector(
    const CableSpan::Impl& cable,
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError)
{
    // Reset path error vector to zero.
    pathError *= 0;

    // The element in the path error vector to write to.
    int pathErrorIx = -1;
    // Iterator over the straight line segments connecting each curve segment.
    auto lineIt = lines.begin();

    cable.forEachActiveCurveSegment(
        s,
        [&](const CurveSegment& curve)
        {
            const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);
            // Compute path error at first contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++pathErrorIx) =
                    calcPathErrorElement(*lineIt, dataPos.X_GP.R(), axis);
            }
            ++lineIt;
            // Compute path error at final contact point.
            for (CoordinateAxis axis : axes) {
                pathError(++pathErrorIx) =
                    calcPathErrorElement(*lineIt, dataPos.X_GQ.R(), axis);
            }
        });
}

} // namespace

//==============================================================================
//                             Path Error Jacobian
//==============================================================================
// This section contains helper functions for computing the jacobian of the
// previously defined path error vector to the natural geodesic corrections.
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
    const Real a                               = dataInst.jacobi_Q[0];
    const Real r                               = dataInst.jacobi_Q[1];

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

template <size_t N>
void calcPathErrorJacobian(
    const CableSpan::Impl& cable,
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J)
{
    // Number of free coordinates for a generic geodesic.
    constexpr int NQ = c_GeodesicDOF;

    const int numberOfCurvesInContact = static_cast<int>(lines.size()) - 1;

    // Reset the values in the jacobian.
    J *= 0.;
    SimTK_ASSERT3(
        J.ncol() == numberOfCurvesInContact * NQ &&
            J.nrow() ==
                numberOfCurvesInContact * MatrixWorkspace::c_NumConstraints,
        "Jacobian matrix size (%d x %d) does not match number of curves (%d)",
        J.nrow(),
        J.ncol(),
        numberOfCurvesInContact);

    // Current indexes to write the elements of the jacobian to,
    int row = 0;
    int col = 0;

    // Helper for writing a computed block to the jacobian.
    auto AddBlock = [&](const Vec4& block, int colOffset = 0)
    {
        for (int ix = 0; ix < 4; ++ix) {
            J.updElt(row, col + colOffset + ix) += block[ix];
        }
    };

    auto linesIt = lines.begin();
    cable.forEachActiveCurveSegment(
        s,
        [&](const CurveSegment& curve)
        {
            const CurveSegmentData::Pos& dataPos = curve.getDataPos(s);
            for (CoordinateAxis axis : axes) {
                const UnitVec3 a_P     = dataPos.X_GP.R().getAxisUnitVec(axis);
                const LineSegment& l_P = *linesIt;

                AddBlock(calcJacobianOfPathErrorAtP(curve, s, l_P, a_P));

                const ObstacleIndex prevObsIx =
                    curve.findPrevObstacleInContactWithCable(s);
                if (prevObsIx.isValid()) {
                    AddBlock(
                        calcJacobianOfNextPathError(
                            cable.getObstacleCurveSegment(prevObsIx),
                            s,
                            l_P,
                            a_P),
                        -NQ);
                }
                ++row;
            }

            ++linesIt;

            for (CoordinateAxis axis : axes) {
                const UnitVec3 a_Q     = dataPos.X_GQ.R().getAxisUnitVec(axis);
                const LineSegment& l_Q = *linesIt;

                AddBlock(calcJacobianOfPathErrorAtQ(curve, s, l_Q, a_Q));

                const ObstacleIndex nextObsIx =
                    curve.findNextObstacleInContactWithCable(s);
                if (nextObsIx.isValid()) {
                    AddBlock(
                        calcJacobianOfPrevPathError(
                            cable.getObstacleCurveSegment(nextObsIx),
                            s,
                            l_Q,
                            a_Q),
                        NQ);
                }
                ++row;
            }

            col += NQ;
        });
}

} // namespace

//==============================================================================
//                      Solving for NaturalGeodesicCorrections
//==============================================================================

namespace
{

// Solve for the geodesic corrections by attempting to set the path error to
// zero. We call this after having filled in the pathError vector and pathError
// jacobian in the MatrixWorkspace. The result is a vector of Corrections for
// each curve.
void calcPathCorrections(MatrixWorkspace& data)
{
    // TODO add explanation...
    // Add a cost to changing the length.
    for (int i = 0; i < data.nObstaclesInContact; ++i) {
        int r = data.nObstaclesInContact *
                    MatrixWorkspace::c_NumPathErrorConstraints + i;
        int c = c_GeodesicDOF * (i + 1) - 1;
        data.pathErrorJacobian.set(r, c, data.maxPathError);
    }

    data.inverse = data.pathErrorJacobian;
    data.inverse.solve(data.pathError, data.pathCorrection);
    data.pathCorrection *= -1.;
}

// Given a correction vector computed from minimizing the cost function,
// compute the maximum allowed stepsize along that correction vector.
// The allowed step is computed by approximating the surface locally as
// a circle in each direction, and limiting the radial displacement on that
// circle.
void calcMaxAllowedCorrectionStepSize(
    const CurveSegment& curve,
    const State& s,
    const NaturalGeodesicCorrection& c,
    Real maxAngularStepsize,
    Real& maxAllowedStepSize)
{
    // Helper to update the maximum allowed stepsize such that the linear
    // displacement does not violate the max allowed step size.
    auto UpdateMaxStepSize = [&](Real maxDisplacementEstimate, Real curvature)
    {
        // Use the local curvature to convert the allowed angular step size to
        // a linear step size.
        const Real maxAllowedDisplacement = maxAngularStepsize / curvature;
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

} // namespace

//------------------------------------------------------------------------------
//                      CableSpan::Impl Cache Computation
//------------------------------------------------------------------------------

const MatrixWorkspace& CableSpan::Impl::calcDataInst(const State& s) const
{
    CableSpanData::Instance& dataInst = updDataInst(s);

    // Axes considered when computing the path error.
    const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

    // Make sure all curve segments are realized to position stage.
    // This will transform all last computed geodesics to Ground frame, and
    // will update each curve's WrappingStatus.
    for (ObstacleIndex ix(0); ix < getNumObstacles(); ++ix) {
        getObstacleCurveSegment(ix).calcDataPos(s);
    }

    // Now that the WrappingStatus of all curve segments is known: Count
    // the number of obstacles in contact with the path.
    const int nActive = countActive(s);

    // If some obstacles are in contact with the cable the path error needs
    // to be checked. If the path error is small, i.e. there are no "kinks"
    // anywhere, the current path is OK. If the path error is too large,
    // corrections need to be computed for each curve segment in order to
    // drive the path error to zero.

    // Grab the shared data cache for helping with computing the path
    // corrections. This data is only used as an intermediate variable, and
    // will be discarded after each iteration. Note that the number active
    // segments determines the sizes of the matrices involved.
    MatrixWorkspace& data = updDataInst(s).updOrInsert(nActive);

    // Reset before computing new values.
    data.maxPathError = 0.;

    // Compute the straight-line segments of this cable span.
    calcLineSegments(
        *this,
        s,
        calcOriginPointInGround(s),
        calcTerminationPointInGround(s),
        data.lineSegments);

    // If the path contains no curved segments it is a straight line.
    if (nActive == 0) {
        return data;
    }

    // Evaluate the current path error as the misalignment of the straight
    // line segments with the curve segment's tangent vectors at the
    // contact points.
    calcPathErrorVector<2>(*this, s, data.lineSegments, axes, data.pathError);
    data.maxPathError = data.pathError.normInf();

    // Only proceed with computing the jacobian and geodesic corrections if the
    // path error is large.
    if (data.maxPathError < getParameters().smoothnessTolerance) {
        data.pathCorrection *= 0.;
        return data;
    }

    // Evaluate the path error jacobian to the natural geodesic corrections
    // of each curve segment.
    calcPathErrorJacobian<2>(
        *this,
        s,
        data.lineSegments,
        axes,
        data.pathErrorJacobian);

    // Compute the geodesic corrections for each curve segment: This gives
    // us a correction vector in a direction that reduces the path error.
    calcPathCorrections(data);

    // Compute the maximum allowed step size that we take along the
    // correction vector.
    Real stepSize     = 1.;
    int activeCurveIx = 0; // Index of curve in all that are in contact.
    forEachActiveCurveSegment(
        s,
        [&](const CurveSegment& curve)
        {
            // Each curve segment will evaluate if the stepsize is not too large
            // given the local curvature.
            calcMaxAllowedCorrectionStepSize(
                curve,
                s,
                data.getCurveCorrection(activeCurveIx),
                getParameters().solverMaxStepSize,
                stepSize);
            ++activeCurveIx;
        });
    // Apply the maximum stepsize to the corrections.
    data.pathCorrection *= stepSize;

    return data;
}

//==============================================================================
//                         CableSubsystemTestHelper
//==============================================================================
// Test correctness of the path error jacobian by applying a small perturbation
// and comparing the effect on the path error. The perturbation is applied in
// the form of a small geodesic correction along each degree of freedom
// subsequently.
// Only curve segments that are in contact with their respective obstacles are
// tested.
bool CableSubsystemTestHelper::applyPerturbationTest(
    const MultibodySystem& system,
    const CableSubsystem& subsystem,
    const State& s,
    Real perturbation,
    Real bound,
    std::ostream& os)
{
    // Result of this perturbation test.
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

        // Apply a small perturbation to each of the curve segments in the form
        // of a geodesic correction. There are 4 DOF for each geodesic, and
        // just to be sure we apply a negative and positive perturbation along
        // each DOF of each geodesic.
        for (int i = 0; i < (c_GeodesicDOF * 2 * nActive); ++i) {
            // We do not want to mess with the actual state, so we make a copy.
            const State sCopy = s;
            system.realize(sCopy, Stage::Position);

            // Trigger realizing position level cache, resetting the
            // configuration.
            const CableSpanData::Pos& dataPos = cable.getDataPos(sCopy);

            // The wrapping status of each obstacle is reset after copying the
            // state. The matrix sizes are no longer compatible if this results
            // in a change in wrapping status (if we are touching down during
            // this realization for example). We can only test the jacobian if
            // the wrapping status of all obstacles remains the same.
            if (cable.countActive(sCopy) != nActive) {
                os << "Cable wrapping status changed after cloning the state: "
                      "Skipping perturbation test\n";
                break;
            }

            MatrixWorkspace& data = cable.updDataInst(s).updOrInsert(nActive);

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

            calcPathErrorVector<2>(
                cable,
                sCopy,
                data.lineSegments,
                axes,
                data.pathError);

            calcPathErrorJacobian<2>(
                cable,
                sCopy,
                data.lineSegments,
                axes,
                data.pathErrorJacobian);

            Vector predictedPathError =
                data.pathErrorJacobian * data.pathCorrection + data.pathError;

            std::vector<ObstacleWrappingStatus> prevWrappingStatus{};
            for (const CurveSegment& curve : cable.m_curveSegments) {
                prevWrappingStatus.push_back(
                    curve.getDataInst(sCopy).wrappingStatus);
            }

            int activeCurveIx = 0; // Index of curve in all that are in contact.
            cable.forEachActiveCurveSegment(
                sCopy,
                [&](const CurveSegment& curve)
                {
                    curve.applyGeodesicCorrection(
                        sCopy,
                        data.getCurveCorrection(activeCurveIx));
                    ++activeCurveIx;
                });

            for (const CurveSegment& curve : cable.m_curveSegments) {
                curve.invalidatePosEntry(sCopy);
            }

            for (ObstacleIndex ix(0); ix < cable.getNumObstacles(); ++ix) {
                cable.getObstacleCurveSegment(ix).calcDataPos(sCopy);
            }

            calcLineSegments(
                cable,
                sCopy,
                dataPos.originPoint_G,
                dataPos.terminationPoint_G,
                data.lineSegments);
            calcPathErrorVector<2>(
                cable,
                sCopy,
                data.lineSegments,
                axes,
                data.pathError);

            // We can only use the path error jacobian if the small
            // perturbation did not trigger any liftoff or touchdown on any
            // obstacles. If any CurveSegment's WrappingStatus has changed, we
            // will not continue with the test. Since the perturbation is
            // small, this is unlikely to happen often.
            bool WrappingStatusChanged = false;
            {
                int ix = -1;
                for (const CurveSegment& curve : cable.m_curveSegments) {
                    WrappingStatusChanged |=
                        prevWrappingStatus.at(++ix) !=
                        curve.getDataInst(sCopy).wrappingStatus;
                }
            }
            if (WrappingStatusChanged) {
                os << "Wrapping status changed: Stopping test\n";
                break;
            }

            const Real predictionError =
                (predictedPathError - data.pathError).norm();

            bool passedTest = std::abs(predictionError / perturbation) <= bound;
            if (!passedTest) {
                os << "FAILED perturbation test for correction = "
                   << data.pathCorrection << "\n";
                os << "    Got      : " << data.pathError << "\n";
                os << "    Predicted: " << predictedPathError << "\n";
                os << "    Err      : "
                   << (predictedPathError - data.pathError) / perturbation
                   << "\n";
                os << "    Err norm : "
                   << (predictedPathError - data.pathError).normInf() /
                          perturbation
                   << "\n";
            } else {
                os << "PASSED perturbation test for correction = "
                   << data.pathCorrection << "\n";
                os << " ( error = "
                   << (predictedPathError - data.pathError).normInf() /
                          perturbation
                   << " )\n";
            }
            success &= passedTest;
        }
    }
    return success;
}

//==============================================================================
//                          CableSubsystem::Impl
//==============================================================================

const CableSpan::Impl& CableSubsystem::Impl::getCableImpl(
    CableSpanIndex index) const
{
    return getCable(index).getImpl();
}

int CableSubsystem::Impl::realizeSubsystemTopologyImpl(State& state) const
{
    // Briefly allow writing into the Topology cache. After this the
    // Topology cache is const.
    Impl* mutableThis = const_cast<Impl*>(this);

    CableSpanData::Instance dataInst{};
    CacheEntryIndex indexDataInst = allocateCacheEntry(
        state,
        Stage::Instance,
        Stage::Infinity,
        new Value<CableSpanData::Instance>(dataInst));

    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        CableSpan& path = mutableThis->updCable(ix);
        path.updImpl().realizeTopology(state, indexDataInst);
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
        cable.getImpl().calcDecorativeGeometryAndAppend(state, decorations);
    }
    return 0;
}

//==============================================================================
//                          CableSubsystem
//==============================================================================

CableSubsystem::CableSubsystem(MultibodySystem& mbs)
{
    adoptSubsystemGuts(new Impl());
    mbs.adoptSubsystem(*this);
}

const MultibodySystem& CableSubsystem::getMultibodySystem() const
{
    return getImpl().getMultibodySystem();
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

const CableSubsystem::Impl& CableSubsystem::getImpl() const
{
    return SimTK_DYNAMIC_CAST_DEBUG<const Impl&>(getSubsystemGuts());
}

CableSubsystem::Impl& CableSubsystem::updImpl()
{
    return SimTK_DYNAMIC_CAST_DEBUG<Impl&>(updSubsystemGuts());
}

//==============================================================================
//                                CableSpan
//==============================================================================

CableSpan::CableSpan() : m_Impl(new Impl())
{}

CableSpan::CableSpan(
    CableSubsystem& subsystem,
    MobilizedBodyIndex originBody,
    const Vec3& originPoint_B,
    MobilizedBodyIndex terminationBody,
    const Vec3& terminationPoint_B) :
    m_Impl(new Impl(
        originBody,
        originPoint_B,
        terminationBody,
        terminationPoint_B))
{
    CableSpanIndex ix = subsystem.updImpl().adoptCable(*this);
    updImpl().setSubsystem(subsystem, ix);
}

ObstacleIndex CableSpan::addObstacle(
    MobilizedBodyIndex obstacleBody,
    const Transform& X_BS,
    std::shared_ptr<const ContactGeometry> obstacleGeometry)
{
    return updImpl().addObstacle(
        obstacleBody,
        X_BS,
        std::move(obstacleGeometry),
        Vec3{NaN}); // Invalid contactPointHint_S
}

ObstacleIndex CableSpan::addObstacle(
    MobilizedBodyIndex obstacleBody,
    const Transform& X_BS,
    std::shared_ptr<const ContactGeometry> obstacleGeometry,
    const Vec3& contactPointHint_S)
{
    return updImpl().addObstacle(
        obstacleBody,
        X_BS,
        std::move(obstacleGeometry),
        contactPointHint_S);
}

int CableSpan::getNumObstacles() const
{
    return getImpl().getNumObstacles();
}

const MobilizedBodyIndex& CableSpan::getObstacleMobilizedBodyIndex(
    ObstacleIndex ix) const
{
    return getImpl().getObstacleCurveSegment(ix).getMobilizedBodyIndex();
}

void CableSpan::setObstacleMobilizedBodyIndex(
    ObstacleIndex ix,
    MobilizedBodyIndex body)
{
    updImpl().updObstacleCurveSegment(ix).setMobilizedBodyIndex(body);
}

const Transform& CableSpan::getObstacleXformSurfaceToBody(
    ObstacleIndex ix) const
{
    return getImpl().getObstacleCurveSegment(ix).getXformSurfaceToBody();
}
void CableSpan::setObstacleXformSurfaceToBody(
    ObstacleIndex ix,
    const Transform& X_BS)
{
    updImpl().updObstacleCurveSegment(ix).setXformSurfaceToBody(X_BS);
}

const ContactGeometry& CableSpan::getObstacleContactGeometry(
    ObstacleIndex ix) const
{
    return getImpl().getObstacleCurveSegment(ix).getContactGeometry();
}
void CableSpan::setObstacleContactGeometry(
    ObstacleIndex ix,
    std::shared_ptr<const ContactGeometry> geometry)
{
    updImpl().updObstacleCurveSegment(ix).setContactGeometry(geometry);
}

Vec3 CableSpan::getObstacleContactPointHint(ObstacleIndex ix) const
{
    return getImpl().getObstacleCurveSegment(ix).getContactPointHint();
}

void CableSpan::setObstacleContactPointHint(
    ObstacleIndex ix,
    Vec3 initialContactPointHint)
{
    updImpl().updObstacleCurveSegment(ix).setContactPointHint(
        initialContactPointHint);
}

bool CableSpan::isInContactWithObstacle(const State& state, ObstacleIndex ix)
    const
{
    getImpl().realizePosition(state);
    return getImpl().getObstacleCurveSegment(ix).isInContactWithSurface(state);
}

Transform CableSpan::calcCurveSegmentInitialFrenetFrame(
    const State& state,
    CableSpanObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl().getObstacleCurveSegment(ix).getDataPos(state).X_GP;
}

Transform CableSpan::calcCurveSegmentFinalFrenetFrame(
    const State& state,
    CableSpanObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl().getObstacleCurveSegment(ix).getDataPos(state).X_GQ;
}

Real CableSpan::calcCurveSegmentArcLength(
    const State& state,
    CableSpanObstacleIndex ix) const
{
    getImpl().realizePosition(state);
    return getImpl().getObstacleCurveSegment(ix).getDataInst(state).arcLength;
}

Real CableSpan::getCurveSegmentAccuracy() const
{
    return getImpl().getParameters().geodesicIntegratorAccuracy;
}

void CableSpan::setCurveSegmentAccuracy(Real accuracy)
{
    updImpl().updParameters().geodesicIntegratorAccuracy = accuracy;
}

int CableSpan::getSolverMaxIterations() const
{
    return getImpl().getParameters().solverMaxIterations;
}

void CableSpan::setSolverMaxIterations(int maxIterations)
{
    updImpl().updParameters().solverMaxIterations = maxIterations;
}

Real CableSpan::getSmoothnessTolerance() const
{
    return getImpl().getParameters().smoothnessTolerance;
}

void CableSpan::setSmoothnessTolerance(Real tolerance)
{
    updImpl().updParameters().smoothnessTolerance = tolerance;
}

Real CableSpan::calcLength(const State& s) const
{
    return getImpl().getDataPos(s).cableLength;
}

Real CableSpan::calcLengthDot(const State& s) const
{
    return getImpl().getDataVel(s).lengthDot;
}

void CableSpan::applyBodyForces(
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    return getImpl().applyBodyForces(s, tension, bodyForcesInG);
}

void CableSpan::calcDecorativePathPoints(
    const State& state,
    const std::function<void(Vec3 point_G)>& sink) const
{
    getImpl().calcDecorativePathPoints(state, sink);
}

Real CableSpan::calcCablePower(const State& state, Real tension) const
{
    return getImpl().calcCablePower(state, tension);
}

int CableSpan::getNumSolverIterations(const State& state) const
{
    return getImpl().getDataPos(state).loopIter;
}

Real CableSpan::getSmoothness(const State& state) const
{
    return getImpl().getDataPos(state).smoothness;
}
