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

#include "CableSpan_SubsystemTestHelper_Impl.h"

#include "simbody/internal/common.h"
#include "simbody/internal/CableSpan.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include "SimTKcommon/internal/State.h"

using namespace SimTK;

// Define a unique index for cable segments.
namespace {
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSegmentIndex);
}

// For convenience.
using ObstacleIndex = CableSpanObstacleIndex;
using ViaPointIndex = CableSpanViaPointIndex;

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
/* The Frenet frame plays an important role in defining the state of a geodesic.
In the end it is simply a Transform, which is why we define it as an alias.

The position of the Frenet frame is defined to lie on the geodesic, but for
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
    // Given the number of CurveSegments that are in contact with their
    // respective obstacle's surface, contruct a MatrixWorkspace of correct
    // dimensions.
    explicit MatrixWorkspace(int problemSize)
    {
        // Total number of obstacles in contact with cable.
        numObstaclesInContact = problemSize;
        // Dimension of the free variable vector: all geodesic corrections.
        const int nQ = numObstaclesInContact * c_GeodesicDOF;
        // Total number of constraints (2 per obstacle).
        const int nC = numObstaclesInContact * 2;

        lineSegments.resize(numObstaclesInContact + 1);
        pathCorrection = Vector(nQ, 0.);

        normalPathError   = Vector(nC, 0.);
        binormalPathError = Vector(nC, 0.);

        normalPathErrorJacobian = Matrix(nC, nQ, 0.);

        A      = Matrix(nC, nC, 0.);
        b      = Vector(nC, 0.);
        lambda = Vector(nC, 0.);

        lengthGradient               = Vector(nQ, 0.);
        lengthHessian                = Matrix(nQ, nQ, 0.);
        lengthHessianInverseEstimate = Matrix(nQ, nQ, 0.);

        singularValues      = Vector(nQ, 0.);
        leftSingularValues  = Matrix(nQ, nQ, 0.);
        rightSingularValues = Matrix(nQ, nQ, 0.);
    }

    // Return the NaturalGeodesicCorrection for the curve segment at the
    // "active" curve index, where "active" means counting those that are in
    // contact with the obstacle.
    static NaturalGeodesicCorrection getCurveCorrection(
        const Vector& pathCorrection,
        int activeCurveIx)
    {
        const int eltIx = activeCurveIx * c_GeodesicDOF;
        SimTK_ASSERT(0 <= eltIx && eltIx + 3 < pathCorrection.size(), "Invalid size of pathCorrection vector.");
        return {
            pathCorrection[eltIx],
            pathCorrection[eltIx + 1],
            pathCorrection[eltIx + 2],
            pathCorrection[eltIx + 3]};
    }

    /* The number of obstacles in contact with the cable. */
    int numObstaclesInContact;
    /* All straight line segments in the current path. */
    std::vector<LineSegment> lineSegments;
    /* The path correction vector contains the NaturalGeodesicCorrection
    vector of each active obstacle stacked as a vector. This vector is
    computed at each solver iteration, and applying these corrections
    attempts to drive the path error vector to zero. */
    Vector pathCorrection;
    /* The normal path error vector is a measure of the angle error between the
    straight line segments and the surface normals at their contact points.

    For each active obstacle two normal path error elements are computed:
    1. dot(e_P, n_P)
    2. dot(e_Q, n_Q)
    where the subscripts P, Q denote the first and final contact frames, e
    denotes the straight line segment's direction vector, and n denotes the
    surface normal.

    Stacking all path error elements, of all active obstacles, gives the
    normal path error vector. */
    Vector normalPathError;
    /* The binormal path error vector is defined analogously to the normal path
    error vector, but measures the angle error against the binormal vector.

    For each active obstacle two normal path error elements are computed:
    1. dot(e_P, b_P)
    2. dot(e_Q, b_Q)
    where b denotes the binormal direction.

    Stacking all path error elements, of all active obstacles, gives the
    binormal path error vector. */
    Vector binormalPathError;
    /* The maximum of the normal and binormal path errors. */
    Real maxPathError = NaN;
    /* The Jacobian of the normal path error vector to the natural geodesic
    corrections of all active obstacles. */
    Matrix normalPathErrorJacobian;

    /* The Lagrangian multipliers used to drive the normal path errors to zero.
    */
    Vector lambda;
    /* The matrix A required for computing the Lagrangian multipliers λ from
    A * λ = b */
    Matrix A;
    /* The vector b required for computing the Lagrangian multipliers λ from
    A * λ = b */
    Vector b;
    /* The factorization used to compute the Lagrangian multipliers. */
    FactorQTZ factor;

    /* The gradient of the total cable length to the natural geodesic
    corrections of each active curve segment. */
    Vector lengthGradient;
    /* The Hessian of the total cable length to the natural geodesic corrections
    of each active curve segment. */
    Matrix lengthHessian;
    /* The inverse of the estimate of the Hessian used during optimization. */
    Matrix lengthHessianInverseEstimate;

    /* The Singular Value Decomposition of the Hessian of the length. */
    FactorSVD svd;
    Vector singularValues;
    Matrix leftSingularValues;
    Matrix rightSingularValues;

    /* This flag indicates if the optimizer stopping criterium was met. If so,
    the pathCorrection vector was not computed, because the optimal path is
    reached. */
    bool converged = false;
};

//------------------------------------------------------------------------------
//  CableSpanData
//------------------------------------------------------------------------------
/* CableSpanData is a data structure used by class CableSpan::Impl to store
relevant quantities in the State's cache.

The struct has three members Instance, Position and Velocity that correspond to
the stages (Stage::Instance, Stage::Position and Stage::Velocity) at which they
can be computed.

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
            for (size_t i = matrixWorkspaces.size();
                 i <= static_cast<size_t>(nActive);
                 ++i) {
                matrixWorkspaces.emplace_back(static_cast<int>(i));
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
    struct Position final {
        // Position of the cable origin point in the ground frame.
        Vec3 originPoint_G{NaN, NaN, NaN};
        // Position of the cable termination point in the ground frame.
        Vec3 terminationPoint_G{NaN, NaN, NaN};
        // Tangent vector to the cable at the origin point, in the ground frame.
        UnitVec3 originTangent_G{NaN, NaN, NaN};
        // Tangent vector to the cable at the termination point, in the ground
        // frame.
        UnitVec3 terminationTangent_G{NaN, NaN, NaN};
        // Tangent vectors to the cable segments incoming to (i.e., immediately
        // prior to) each via point in the path, in the ground frame.
        Array_<UnitVec3, ViaPointIndex> viaPointInTangents_G;
        // Tangent vectors to the cable segments outgoing from (i.e.,
        // immediately following) each via point in the path, in the ground frame.
        Array_<UnitVec3, ViaPointIndex> viaPointOutTangents_G;
        // The total cable length.
        Real cableLength = NaN;
        // The path smoothness, computed as the maximum angular misalignment at
        // the points where the straight- and curved-line segments meet,
        // measured in radians.
        Real smoothness = NaN;
        // Number of iterations the solver used to compute the cable's path.
        int loopIter = -1;
    };
    struct Velocity final {
        // Time derivative of the total cable length.
        Real lengthDot = NaN;
    };
};

//------------------------------------------------------------------------------
//  CurveSegmentData
//------------------------------------------------------------------------------
/* CurveSegmentData is a data structure used by class CurveSegment to store
relevent quantities in the State's cache.

The struct has two member structs Instance and Position that correspond to the
stages (Stage::Instance and Stage::Position) at which they can be computed.

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
        std::vector<ContactGeometry::GeodesicKnotPoint> geodesicKnotPoints;
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
    struct Position final {
        // Position and orientation of contact geometry w.r.t. ground.
        Transform X_GS;
        // Frenet frame at the initial contact point w.r.t. ground.
        FrenetFrame X_GP;
        // Frenet frame at the final contact point w.r.t. ground.
        FrenetFrame X_GQ;
        // Variation of the Frenet frame position at frame P. Each column is
        // the variation of the position to each of the natural geodesic
        // variations.
        //
        // For example, let x_P be the position of the Frenet frame P, and q
        // the geodesic variation, then:
        //
        // dx_P / dq = v_P
        Mat34 v_P;
        // Variation of the Frenet frame orientation at frame P. Each column is
        // the variation of the orientation to each of the natural geodesic
        // variations, as a rotation vector.
        //
        // Note that the columns of v_P and w_P together are simply the spatial
        // vectors describing the variation of the Frenet frame P.
        //
        // For example, let S(x) be the 3x3 skew-symmetric matrix related to
        // the cross product by x % y = S(x) * y, then
        //
        // dt_P / dq = -S(t_P) * w_P
        // dn_P / dq = -S(n_P) * w_P
        // db_P / dq = -S(b_P) * w_P
        //
        // where t_P, n_P, b_P denote the tangent, normal and binormal
        // directions.
        Mat34 w_P;
        // Variation of the Frenet frame position at frame Q.
        Mat34 v_Q;
        // Variation of the Frenet frame orientation at frame Q.
        Mat34 w_Q;
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
when computing a new geodesic over the obstacle's surface. */
struct IntegratorTolerances {
    // The accuracy used to control the stepsize of the variable step geodesic
    // integrator.
    Real geodesicIntegratorAccuracy = 1e-9;
    // The tolerance used when projecting the geodesic's state to the surface
    // during integration.
    // TODO the surface projection tolerance should be linked to the geodesic
    // integrator accuracy.
    Real constraintProjectionTolerance = 1e-9;
    // TODO this is not connected to anything yet.
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
    // The algorithm used to compute the optimal path.
    CableSpanAlgorithm algorithm = CableSpanAlgorithm::MinimumLength;
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

// Compute the Frenet frame at a knot point of a geodesic.
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

// Evaluate whether a point lies below the given surface.
bool isPointBelowSurface(const ContactGeometry& geometry, const Vec3& point_S)
{
    // We are not below the surface if the surface is not defined at the point.
    if (!geometry.isSurfaceDefined(point_S)) {
        return false;
    }
    // NOTE surface value is negative above surface.
    return geometry.calcSurfaceValue(point_S) > 0.;
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

    int realizeSubsystemModelImpl(State& state) const override;

    int realizeSubsystemInstanceImpl(const State& state) const override;

    int realizeSubsystemTimeImpl(const State& state) const override;

    int realizeSubsystemPositionImpl(const State& state) const override;

    int realizeSubsystemVelocityImpl(const State& state) const override;

    int realizeSubsystemDynamicsImpl(const State& state) const override;

    int realizeSubsystemAccelerationImpl(const State& state) const override;

    int realizeSubsystemReportImpl(const State& state) const override;

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
//                               Curve Segment
//==============================================================================
/* The CableSpan's path consists of straight line segments between obstacles and
via points, and curved segments over the obstacles. This struct represents the
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
        CableSegmentIndex cableSegmentIndex,
        ObstacleIndex obstacleIndex,
        MobilizedBodyIndex obstacleBody,
        Transform X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry,
        const Vec3& contactPointHint_S) :
        m_subsystem(subsystem),
        m_cableIndex(cableIndex),
        m_cableSegmentIndex(cableSegmentIndex),
        m_obstacleIndex(obstacleIndex),
        m_body(obstacleBody),
        m_X_BS(std::move(X_BS)),
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
            dataInst.wrappingStatus = ObstacleWrappingStatus::LiftedFromSurface;
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
            new Value<CurveSegmentData::Position>());
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

    const CurveSegmentData::Position& getDataPos(const State& state) const
    {
        return Value<CurveSegmentData::Position>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataPos));
    }

    CurveSegmentData::Position& updDataPos(const State& state) const
    {
        return Value<CurveSegmentData::Position>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataPos));
    }

    //--------------------------------------------------------------------------
    // Model configuration accessors.
    //--------------------------------------------------------------------------

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

    const Transform& getTransformSurfaceToBody() const
    {
        return m_X_BS;
    }

    void setTransformSurfaceToBody(const Transform& X_BS)
    {
        m_X_BS = X_BS;
    }

    const CableSubsystem& getSubsystem() const
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CurveSegment::getSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

        return *m_subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CurveSegment::updSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

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

    CableSegmentIndex getCableSegmentIndex() const
    {
        return m_cableSegmentIndex;
    }

    const IntegratorTolerances& getIntegratorTolerances() const;

    //--------------------------------------------------------------------------
    // Utility functions.
    //--------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------
    //  Helper Functions: Finding next and previous CurveSegments
    //--------------------------------------------------------------------------
    // The functions below essentially filter any obstacles that are not in
    // contact with the cable.

    // Find the first obstacle before this segment that is in contact with the
    // cable. Returns an invalid index if there is none, e.g., the previous
    // point in the path is a via point or the path origin.
    ObstacleIndex findPrevObstacleInContactWithCable(const State& state) const;
    // Find the first obstacle after this segment that is in contact with the
    // cable. Returns an invalid index if there is none, e.g., the next point
    // in the path is a via point or the path termination.
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
            !isPointBelowSurface(getContactGeometry(), prevPoint_S),
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
            !isPointBelowSurface(getContactGeometry(), nextPoint_S),
            "CurveSegment::Impl::assertSurfaceBounds",
            "Next point lies inside the surface");
        return nextPoint_S;
    }

    //--------------------------------------------------------------------------
    //  Updating CurveSegmentData::Instance Helpers
    //--------------------------------------------------------------------------

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
            // For analytic surfaces we can compute the inital and final Frenet
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
            dataInst.geodesicKnotPoints.clear();
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
                    dataInst.geodesicKnotPoints.push_back(q);
                });

            // Copy the state at the boundary frames.
            SimTK_ASSERT(
                dataInst.geodesicKnotPoints.size() > 0,
                "CurveSegment: Failed to shoot a new geodesic: buffer empty");
            geodesicBoundaryStates = {
                dataInst.geodesicKnotPoints.front(),
                dataInst.geodesicKnotPoints.back(),
            };
        }

        // Use the geodesic state at the boundary frames to compute the fields
        // in CurveSegmentData::Instance.

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
                static_cast<int>(dataInst.geodesicKnotPoints.size()) - 1;

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
                    dataInst.geodesicKnotPoints.at(1).arcLength -
                    dataInst.geodesicKnotPoints.front().arcLength;
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
                        dataInst.geodesicKnotPoints.at(2).arcLength -
                        dataInst.geodesicKnotPoints.at(1).arcLength;
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
        bool touchdownDetected =
            getContactGeometry().calcNearestPointOnLineImplicitly(
                prevPoint_S,
                nextPoint_S,
                tols.constraintProjectionMaxIterations,
                tols.constraintProjectionTolerance,
                pointOnLineNearSurface_S) ==
            ContactGeometry::NearestPointOnLineResult::PointFallsBelowSurface;

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

    //--------------------------------------------------------------------------
    //  Computing cached data
    //--------------------------------------------------------------------------

    // Apply the correction to the initial condition of the geodesic, and
    // shoot a new geodesic, updating the cache variable.
    void applyGeodesicCorrection(
        const State& state,
        const NaturalGeodesicCorrection& correction) const
    {
        // Copy the correction vector.
        NaturalGeodesicCorrection c = correction;
        // Get the current geodesic data.
        const CurveSegmentData::Instance& dataInst = getDataInst(state);

        // Take the length correction, and add to the current length.
        const Real dl = c[3]; // Length increment is the last element.
        const Real l  = dataInst.arcLength;
        // Clamp the length to be nonnegative.
        const bool clamped         = l + dl < 0.;
        const Real correctedLength = clamped ? 0. : l + dl;
        if (clamped) {
            // Equally distribute the clamping effect over the P and Q frame.
            c[0] += 0.5 * (l + dl);
        }

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

    const CurveSegmentData::Position& calcDataPos(const State& state) const
    {
        // Make sure the Instance level data is up to date.
        const CurveSegmentData::Instance& dataInst = calcDataInst(state);

        // Update the Stage::Position level cache.
        CurveSegmentData::Position& dataPos = updDataPos(state);

        // Store transform from local surface frame to ground.
        dataPos.X_GS = calcSurfaceFrameInGround(state);

        // Check if the obstacle is in contact with the cable.
        if (isInContactWithSurface(state)) {
            // Store the geodesic's Frenet frames in ground frame.
            dataPos.X_GP = dataPos.X_GS.compose(dataInst.X_SP);
            dataPos.X_GQ = dataPos.X_GS.compose(dataInst.X_SQ);

            // Compute the initial Frenet frame variation.
            {
                const FrenetFrame& X_GP = dataPos.X_GP;

                const Real tau = dataInst.torsion_P;
                const Real kn  = dataInst.curvatures_P[0];
                const Real kb  = dataInst.curvatures_P[1];

                dataPos.v_P.col(0) = Vec3(getTangent(X_GP));
                dataPos.v_P.col(1) = Vec3(getBinormal(X_GP));
                dataPos.v_P.col(2) = Vec3(0.);
                dataPos.v_P.col(3) = Vec3(0.);

                dataPos.w_P.col(0) = X_GP.R() * Vec3(tau, 0., kn);
                dataPos.w_P.col(1) = X_GP.R() * Vec3(-kb, 0., -tau);
                dataPos.w_P.col(2) = X_GP.R() * Vec3(0., -1., 0.);
                dataPos.w_P.col(3) = Vec3(0);
            }

            // Compute the final Frenet frame variation.
            {
                const FrenetFrame& X_GQ = dataPos.X_GQ;

                const Real tau = dataInst.torsion_Q;
                const Real kn  = dataInst.curvatures_Q[0];
                const Real kb  = dataInst.curvatures_Q[1];

                const Real a    = dataInst.jacobi_Q[0];
                const Real aDot = dataInst.jacobiDot_Q[0];

                const Real r    = dataInst.jacobi_Q[1];
                const Real rDot = dataInst.jacobiDot_Q[1];

                dataPos.v_Q.col(0) = Vec3(getTangent(X_GQ));
                dataPos.v_Q.col(1) = a * getBinormal(X_GQ);
                dataPos.v_Q.col(2) = r * getBinormal(X_GQ);
                dataPos.v_Q.col(3) = Vec3(getTangent(X_GQ));

                dataPos.w_Q.col(0) = X_GQ.R() * Vec3(tau, 0., kn);
                dataPos.w_Q.col(1) = X_GQ.R() * Vec3(-a * kb, -aDot, -a * tau);
                dataPos.w_Q.col(2) = X_GQ.R() * Vec3(-r * kb, -rDot, -r * tau);
                dataPos.w_Q.col(3) = dataPos.w_Q.col(0);
            }
        } else {
            // No contact = no geodesic.
            dataPos.X_GP.setToNaN();
            dataPos.X_GQ.setToNaN();

            dataPos.v_P.setToNaN();
            dataPos.w_P.setToNaN();
            dataPos.v_Q.setToNaN();
            dataPos.w_Q.setToNaN();
        }

        getSubsystem().markCacheValueRealized(state, m_indexDataPos);

        return dataPos;
    }

    //--------------------------------------------------------------------------
    //  Generating Geodesic Points For Visualization
    //--------------------------------------------------------------------------

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
            // NOTE Depending on what is actually desired by the end-user, this
            // might need revision.

            // For analytic surfaces there is no numerical integrator to give a
            // natural interspacing of knots. So instead I just assume you
            // would like some OK granularity. This can be done by interspacing
            // the knots by a certain "angle" (e.g. 5 degrees is pretty fine),
            // and using the radius of curvature to convert the angle to a
            // spatial interspacing of the samples.
            constexpr Real c_AngularSpacing = 0.087; // ~ 5 degrees

            // To convert the angular spacing to a spatial spacing we need the
            // curvature. The curvature is generally not constant along a
            // geodesic, but there are not a lot of analytic surfaces at the
            // moment (sphere and cylinder?). For these we can just look at the
            // curvature at the boundary frames to know the curvature along the
            // entire geodesic.
            const Real maxCurvature_P = max(abs(dataInst.curvatures_P));
            const Real maxCurvature_Q = max(abs(dataInst.curvatures_Q));
            const Real minRadiusOfCurvature =
                1. / std::max(maxCurvature_P, maxCurvature_Q);

            // Estimate the spatial interspacing from the angular spacing and
            // the radius of curvature.
            const Real lengthInterspacing =
                minRadiusOfCurvature * c_AngularSpacing;

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
             dataInst.geodesicKnotPoints) {
            calcPathPointAndWriteToSink(q);
        }
    }

    /* Resample the curve segment at equal length intervals using n=numSamples
    points, and write them to the output sink. */
    void calcResampledPoints(
        const State& state,
        int numSamples,
        const std::function<void(Vec3 point_G)>& sink) const;

private:
    //--------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------

    // Subsystem info.
    CableSubsystem* m_subsystem = nullptr;
    // The index of this CurveSegment's cable in the CableSubsystem.
    CableSpanIndex m_cableIndex = CableSpanIndex::Invalid();
    // The index of the CableSegment that this CurveSegment belongs to.
    CableSegmentIndex m_cableSegmentIndex = CableSegmentIndex::Invalid();
    // The index of this CurveSegment's obstacle in the cable.
    ObstacleIndex m_obstacleIndex = ObstacleIndex::Invalid();

    // MobilizedBody that the curve segment is attached to.
    MobilizedBodyIndex m_body = MobilizedBodyIndex::Invalid();
    // Surface to body transform.
    Transform m_X_BS;

    // Obstacle surface.
    std::shared_ptr<const ContactGeometry> m_geometry;

    // Topology cache.
    CacheEntryIndex m_indexDataPos        = CacheEntryIndex::Invalid();
    DiscreteVariableIndex m_indexDataInst = DiscreteVariableIndex::Invalid();

    // Initial contact point hint used to setup the initial path.
    Vec3 m_contactPointHint_S{NaN, NaN, NaN};

    //--------------------------------------------------------------------------
    // Helper class for unit tests.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                                 Via Point
//==============================================================================
/* This struct represents a zero-length cable segment associated with a via
point in the path. It manages the cached data such as the position of the via
point in the ground frame and the direction of the incoming and outgoing line
segments. */
class ViaPoint final {
public:
    //--------------------------------------------------------------------------
    // Constructors.
    //--------------------------------------------------------------------------

    ViaPoint() = default;

    ViaPoint(
        CableSubsystem* subsystem,
        CableSpanIndex cableIndex,
        ViaPointIndex viaPointIndex,
        MobilizedBodyIndex body,
        const Vec3& station_B) :
        m_subsystem(subsystem),
        m_cableIndex(cableIndex),
        m_viaPointIndex(viaPointIndex),
        m_body(body),
        m_station_B(station_B)
    {}

    //--------------------------------------------------------------------------
    // Realizing and cache access.
    //--------------------------------------------------------------------------

    // Allocate the cache entries.
    void realizeTopology(State& state)
    {
        const Vec3 initVec3{0.};
        const UnitVec3 initUnitVec3; // internally initializes to all-NaN.
        m_indexStation_G = updSubsystem().allocateCacheEntry(
            state,
            Stage::Position,
            Stage::Infinity,
            new Value<Vec3>(initVec3));
    }

    const Vec3& getStation_G(const State& state) const
    {
        return Value<Vec3>::downcast(
            getSubsystem().getCacheEntry(state, m_indexStation_G));
    }

    Vec3& updStation_G(const State& state) const
    {
        return Value<Vec3>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexStation_G));
    }

    //--------------------------------------------------------------------------
    // Model configuration accessors.
    //--------------------------------------------------------------------------

    ViaPointIndex getIndex() const
    {
        return m_viaPointIndex;
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

    const CableSubsystem& getSubsystem() const
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "ViaPoint::getSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

        return *m_subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "ViaPoint::updSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

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

    const Vec3& getStation_B() const
    {
        return m_station_B;
    }

    void setStation_B(const Vec3& station_B)
    {
        m_station_B = station_B;
    }

    //--------------------------------------------------------------------------
    // Computing cached data
    //--------------------------------------------------------------------------

    const Vec3& calcStation_G(const State& state) const
    {
        Vec3& station_G = updStation_G(state);
        station_G = getMobilizedBody().getBodyTransform(
            state).shiftFrameStationToBase(m_station_B);
        getSubsystem().markCacheValueRealized(state, m_indexStation_G);
        return station_G;
    }

    //--------------------------------------------------------------------------
    // Helper function: via point visualization
    //--------------------------------------------------------------------------

    void calcDecorativePathPoints(
        const State& state,
        const std::function<void(Vec3 point_G)>& sink) const
    {
        sink(getStation_G(state));
    }

private:
    //--------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------
    // Subsystem info.
    CableSubsystem* m_subsystem = nullptr;
    // The index of this ViaPoint's cable in the CableSubsystem.
    CableSpanIndex m_cableIndex = CableSpanIndex::Invalid();
    // The index of this ViaPoint in the CableSpan.
    ViaPointIndex m_viaPointIndex = ViaPointIndex::Invalid();

    // MobilizedBody that the via point is attached to.
    MobilizedBodyIndex m_body = MobilizedBodyIndex::Invalid();
    // The position of the via point in the body frame.
    Vec3 m_station_B{NaN};

    // Topology cache.
    CacheEntryIndex m_indexStation_G = CacheEntryIndex::Invalid();
};

//==============================================================================
//                              Cable Segment
//==============================================================================
/* This is a helper class for managing segments of the cable path that are
separated by via points. It stores an ordered list of the obstacle indexes
associated with this segment and indexes to the initial and final via points, if
they exist. An invalid initial via point index indicates that the segment begins
with the cable origin point. Similarly, an invalid final via point index
indicates that the segment ends with the cable termination point. */
class CableSegment {
public:
    CableSegment() = default;

    CableSegment(
        CableSubsystem* subsystem,
        CableSpanIndex cableIndex,
        CableSegmentIndex index) :
        m_subsystem(subsystem),
        m_cableIndex(cableIndex),
        m_index(index)
    {}

    const CableSubsystem& getSubsystem() const
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CableSegment::getSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

        return *m_subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        SimTK_ERRCHK_ALWAYS(
            m_subsystem,
            "CableSegment::updSubsystem()",
            "CableSpan not yet adopted by any CableSubsystem");

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

    CableSegmentIndex getIndex() const
    {
        return m_index;
    }

    void addObstacleIndex(ObstacleIndex obstacleIndex)
    {
        m_obstacleIndexes.push_back(obstacleIndex);
    }

    const Array_<ObstacleIndex>& getObstacleIndexes() const
    {
        return m_obstacleIndexes;
    }

    void setInitialViaPointIndex(ViaPointIndex viaPointIndex)
    {
        m_initialViaPointIndex = viaPointIndex;
    }

    ViaPointIndex getInitialViaPointIndex() const
    {
        return m_initialViaPointIndex;
    }

    void setFinalViaPointIndex(ViaPointIndex viaPointIndex)
    {
        m_finalViaPointIndex = viaPointIndex;
    }

    ViaPointIndex getFinalViaPointIndex() const
    {
        return m_finalViaPointIndex;
    }

private:
    // Subsystem info.
    CableSubsystem* m_subsystem = nullptr;
    // The index of this CableSegment's cable in the CableSubsystem.
    CableSpanIndex m_cableIndex = CableSpanIndex::Invalid();
    // The index of this CableSegment in the CableSpan.
    CableSegmentIndex m_index = CableSegmentIndex::Invalid();

    // The indexes to the wrap obstacles belonging to this segment.
    Array_<ObstacleIndex> m_obstacleIndexes;

    // The indexes to the initial and final via points belonging to this
    // segment.
    ViaPointIndex m_initialViaPointIndex = ViaPointIndex::Invalid();
    ViaPointIndex m_finalViaPointIndex = ViaPointIndex::Invalid();

    //--------------------------------------------------------------------------
    // Helper class for unit tests.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                         Class Cablespan::Impl
//==============================================================================
/* This is the internal implementation class for CableSpan.

The CableSpan::Impl's main job is to set up the solver for computing the path
over the surfaces. A list of CurveSegments is stored, and whenever an obstacle
is added to the CableSpan, this will create a new CurveSegment internally. The
CurveSegment deals with computing the local geodesic, and detecting contact or
liftoff, and manages the related cached data. We then collect info on the
computed geodesics, to evaluate the path smoothness (see
CableSpanData::Instance::pathError). Depending on the given tolerance, we
compute the NaturalGeodesicCorrection vector for each CurveSegment to improve
the smoothness, and try again. */
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
        Vec3 originStation,
        MobilizedBodyIndex terminationBody,
        Vec3 terminationStation) :
        m_originBody(originBody),
        m_originStation(originStation),
        m_terminationBody(terminationBody),
        m_terminationStation(terminationStation)
    {}

    //--------------------------------------------------------------------------
    // Realizing and accessing cache.
    //--------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    // The cached CableSpanData::Instance is shared among those CableSpan
    // managed by one CableSubsystem, and accessed by the provided
    // indexOfDataInstEntry.
    void realizeTopology(State& state, CacheEntryIndex indexOfDataInstEntry)
    {
        m_indexDataInst = indexOfDataInstEntry;

        // Initialize the length of the via point tangent vectors to the number
        // of via points in the cable.
        CableSpanData::Position dataPos;
        const int nViaPoints = getNumViaPoints();
        dataPos.viaPointInTangents_G.resize(nViaPoints);
        dataPos.viaPointOutTangents_G.resize(nViaPoints);
        m_indexDataPos = updSubsystem().allocateCacheEntry(
            state,
            Stage::Position,
            Stage::Infinity,
            new Value<CableSpanData::Position>(dataPos));

        m_indexDataVel = updSubsystem().allocateCacheEntry(
            state,
            Stage::Velocity,
            Stage::Infinity,
            new Value<CableSpanData::Velocity>());

        for (CurveSegment& segment : m_curveSegments) {
            segment.realizeTopology(state);
        }

        for (ViaPoint& viaPoint : m_viaPoints) {
            viaPoint.realizeTopology(state);
        }
    }

    void invalidateTopology()
    {
        if (m_subsystem) {
            getSubsystem().invalidateSubsystemTopologyCache();
        }
    }

    void realizeModel(State& state) const
    {
        // No choices at the moment.
    }

    void realizeInstance(const State& state) const
    {
        // Nothing to compute here yet.
    }

    void realizeTime(const State& state) const
    {
        // Nothing to compute here yet.
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

    void realizeDynamics(const State& state) const
    {
        // Nothing to compute here yet.
    }

    void realizeAcceleration(const State& state) const
    {
        // Nothing to compute here yet.
    }

    void realizeReport(const State& state) const
    {
        // Nothing to compute here yet.
    }

    const CableSpanData::Position& getDataPos(const State& state) const
    {
        realizePosition(state);
        return Value<CableSpanData::Position>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataPos));
    }

    const CableSpanData::Velocity& getDataVel(const State& state) const
    {
        realizeVelocity(state);
        return Value<CableSpanData::Velocity>::downcast(
            getSubsystem().getCacheEntry(state, m_indexDataVel));
    }

    //--------------------------------------------------------------------------
    // Accessors
    //--------------------------------------------------------------------------

    const CurveSegment& getObstacleCurveSegment(ObstacleIndex ix) const
    {
        return m_curveSegments.at(ix);
    }

    CurveSegment& updObstacleCurveSegment(ObstacleIndex ix)
    {
        return m_curveSegments.at(ix);
    }

    const ViaPoint& getViaPoint(ViaPointIndex ix) const
    {
        return m_viaPoints.at(ix);
    }

    ViaPoint& updViaPoint(ViaPointIndex ix)
    {
        return m_viaPoints.at(ix);
    }

    const CableSegment& getCableSegment(CableSegmentIndex ix) const
    {
        return m_cableSegments.at(ix);
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
            m_originStation);
    }

    Vec3 calcTerminationPointInGround(const State& state) const
    {
        return getTerminationBody()
            .getBodyTransform(state)
            .shiftFrameStationToBase(m_terminationStation);
    }

    //--------------------------------------------------------------------------
    // Helper function: Counting active curve segments.
    //--------------------------------------------------------------------------

    // Count the number of CurveSegments that are in contact with the obstacle's
    // surface for a given cable segment.
    int countActive(const State& s, const CableSegment& cableSegment) const
    {
        int counter = 0;
        for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
            const CurveSegment& curve = getObstacleCurveSegment(ix);
            if (curve.isInContactWithSurface(s)) {
                ++counter;
            }
        }
        return counter;
    }

    //--------------------------------------------------------------------------
    // Helper functions: CableSegment initial and final points.
    //--------------------------------------------------------------------------

    Vec3 calcInitialCableSegmentPointInGround(const State& s,
        const CableSegment& cableSegment) const
    {
        if (cableSegment.getInitialViaPointIndex().isValid()) {
            return getViaPoint(cableSegment.getInitialViaPointIndex())
                                .updStation_G(s);
        }
        return calcOriginPointInGround(s);
    }

    Vec3 calcFinalCableSegmentPointInGround(const State& s,
        const CableSegment& cableSegment) const
    {
        if (cableSegment.getFinalViaPointIndex().isValid()) {
            return getViaPoint(cableSegment.getFinalViaPointIndex())
                                .updStation_G(s);
        }
        return calcTerminationPointInGround(s);
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

        SimTK_ASSERT(
            !m_cableSegments.empty(),
            "No cable segments found. An initial cable segment should have "
            "been added by the CableSpan constructor.");

        ObstacleIndex obstacleIx(m_curveSegments.size());

        m_curveSegments.push_back(CurveSegment(
            m_subsystem,
            getIndex(),
            m_cableSegments.back().getIndex(),
            obstacleIx,
            obstacleBody,
            X_BS,
            std::move(obstacleGeometry),
            contactPointHint_S));

        // Add the obstacle to the current CableSegment.
        m_cableSegments.back().addObstacleIndex(obstacleIx);

        return obstacleIx;
    }

    ViaPointIndex addViaPoint(
        MobilizedBodyIndex viaPointBody,
        const Vec3& station_B)
    {
        invalidateTopology();

        SimTK_ASSERT(
            !m_cableSegments.empty(),
            "No cable segments found. An initial cable segment should have "
            "been added by the CableSpan constructor.");

        ViaPointIndex viaPointIx(m_viaPoints.size());

        m_viaPoints.push_back(ViaPoint(
            m_subsystem,
            getIndex(),
            viaPointIx,
            viaPointBody,
            station_B));

        // Mark the end of the last CableSegment with this ViaPoint.
        m_cableSegments.back().setFinalViaPointIndex(viaPointIx);

        // Add a new CableSegment.
        CableSegmentIndex segmentIx(m_cableSegments.size());
        m_cableSegments.push_back(CableSegment(
            m_subsystem,
            getIndex(),
            segmentIx));
        m_cableSegments.back().setInitialViaPointIndex(viaPointIx);

        return viaPointIx;
    }

    MobilizedBodyIndex getOriginBodyIndex() const
    {
        return m_originBody;
    }

    void setOriginBodyIndex(MobilizedBodyIndex originBody)
    {
        m_originBody = originBody;
    }

    MobilizedBodyIndex getTerminationBodyIndex() const
    {
        return m_terminationBody;
    }

    void setTerminationBodyIndex(MobilizedBodyIndex terminationBody)
    {
        m_terminationBody = terminationBody;
    }

    const Vec3& getOriginStation() const
    {
        return m_originStation;
    }

    void setOriginStation(const Vec3& originStation)
    {
        m_originStation = originStation;
    }

    const Vec3& getTerminationStation() const
    {
        return m_terminationStation;
    }

    void setTerminationStation(const Vec3& terminationStation)
    {
        m_terminationStation = terminationStation;
    }

    int getNumObstacles() const
    {
        return m_curveSegments.size();
    }

    int getNumViaPoints() const
    {
        return m_viaPoints.size();
    }

    int getNumCableSegments() const
    {
        return m_cableSegments.size();
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
        const CableSpanData::Position& dataPos = getDataPos(state);
        sink(dataPos.originPoint_G);

        // Write the path points for each cable segment.
        for (CableSegmentIndex ix(0); ix < getNumCableSegments(); ++ix) {
            const CableSegment& cableSegment = getCableSegment(ix);

            // Write the path points for each obstacle in the cable segment.
            for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                const CurveSegment& curve = getObstacleCurveSegment(ix);
                curve.calcDecorativePathPoints(state, sink);
            }

            // Write the path points for the final via point in the cable
            // segment.
            if (cableSegment.getFinalViaPointIndex().isValid()) {
                const ViaPoint& viaPoint = getViaPoint(
                    cableSegment.getFinalViaPointIndex());
                viaPoint.calcDecorativePathPoints(state, sink);
            }
        }

        // Write the path's termination point.
        sink(dataPos.terminationPoint_G);
    }

    void calcCurveSegmentResampledPoints(
        const State& state,
        CableSpanObstacleIndex ix,
        int numSamples,
        const std::function<void(Vec3 point_G)>& sink) const;

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
                    Vec3 surfacePoint_S(NaN);
                    try {
                        surfacePoint_S =
                            curve.getContactGeometry()
                                .projectDownhillToNearestPoint(trackingPoint_S);
                    } catch (...) {
                        // Yes catching all is bad, but since this plotting
                        // option is just for debugging, and this
                        // projection is likely to fail as we get further
                        // from the surface, we just dont plot it if we
                        // don't have it.
                        continue;
                    }
                    const Vec3 surfacePoint_G =
                        X_GS.shiftFrameStationToBase(surfacePoint_S);
                    decorations.push_back(
                        DecorativeLine(trackingPoint_G, surfacePoint_G)
                            .setColor(Orange)
                            .setLineThickness(c_LineThickness));
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
        for (ViaPoint& viaPoint : m_viaPoints) {
            viaPoint.setSubsystem(subsystem);
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

    // Mutable CableSpanData::Position cache access.
    CableSpanData::Position& updDataPos(const State& state) const
    {
        return Value<CableSpanData::Position>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataPos));
    }

    // Mutable CableSpanData::Velocity cache access.
    CableSpanData::Velocity& updDataVel(const State& state) const
    {
        return Value<CableSpanData::Velocity>::updDowncast(
            getSubsystem().updCacheEntry(state, m_indexDataVel));
    }

    //--------------------------------------------------------------------------
    //  Helper function: Creating the first CableSegment
    //--------------------------------------------------------------------------

    void createFirstCableSegment()
    {
        m_cableSegments.push_back(CableSegment(
            m_subsystem,
            getIndex(),
            CableSegmentIndex(0)));
    }

    //--------------------------------------------------------------------------
    //  Computing cached data
    //--------------------------------------------------------------------------

    /* Computes a single correction step for driving the path errors to zero.
    This correction can be computed using different algorithms. The resulting
    path corrections, and intermediate results are written to the
    MatrixWorkspace data output. */
    void calcSolverStep(
        const State& s,
        const CableSegment& cableSegment,
        CableSpanAlgorithm algorithm,
        MatrixWorkspace& data) const;

    const CableSpanData::Position& calcDataPos(const State& s) const;

    void calcDataVel(const State& s, CableSpanData::Velocity& dataVel) const
    {
        const CableSpanData::Position& dataPos = getDataPos(s);

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

        // Compute the lengthening speed of the total cable.
        // Consider a point at x_GQ moving with velocity v_GQ, and similarly a
        // point x_GP moving with velocity v_GP. The line starting at x_GQ and
        // ending at x_GP will have a lenghtening speed of:
        // lengthDot = dot( UnitVec(x_GP - x_GQ), v_GP - v_GQ)
        // The total cable lengthening speed is computed by summing the
        // lengthening speed of all straight line segments (see Scholz2015).
        Real& lengthDot = (dataVel.lengthDot = 0.);

        // Starting point and velocity of straight line segment (the cable
        // origin).
        Vec3 x_GQ = dataPos.originPoint_G;
        Vec3 v_GQ =
            getOriginBody().findStationVelocityInGround(s, m_originStation);

        // Compute the lengthening speed contributed by each cable segment.
        for (CableSegmentIndex ix(0); ix < getNumCableSegments(); ++ix) {
            const CableSegment& cableSegment = getCableSegment(ix);

            // Add the lengthening speed contributed by each obstacle in the
            // cable segment.
            for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                const CurveSegment& curve = getObstacleCurveSegment(ix);
                if (curve.isInContactWithSurface(s)) {
                    const MobilizedBody& mobod = curve.getMobilizedBody();
                    const CurveSegmentData::Position& curveDataPos =
                        curve.getDataPos(s);

                    // Ending point and velocity of straight line segment.
                    const Vec3 x_GP = curveDataPos.X_GP.p();
                    const Vec3 v_GP = CalcPointVelocityInGround(mobod, x_GP);

                    // Compute the lengthening.
                    lengthDot += dot(UnitVec3(x_GP - x_GQ), v_GP - v_GQ);

                    // Set the next straight line starting point and velocity.
                    x_GQ = curveDataPos.X_GQ.p();
                    v_GQ = CalcPointVelocityInGround(mobod, x_GQ);
                }
            }

            // Add the lengthening speed contributed by each via point in the
            // cable segment.
            if (cableSegment.getFinalViaPointIndex().isValid()) {
                const ViaPoint& viaPoint = getViaPoint(
                    cableSegment.getFinalViaPointIndex());
                const MobilizedBody& mobod = viaPoint.getMobilizedBody();

                // Ending point and velocity of straight line segment.
                const Vec3& x_GP = viaPoint.getStation_G(s);
                const Vec3 v_GP = CalcPointVelocityInGround(mobod, x_GP);

                // Compute the lengthening.
                lengthDot += dot(UnitVec3(x_GP - x_GQ), v_GP - v_GQ);

                // Set the next straight line starting point and velocity.
                x_GQ = x_GP;
                v_GQ = v_GP;
            }
        }

        // Ending point and velocity of straight line segment (the cable
        // termination).
        const Vec3 v_GP = getTerminationBody().findStationVelocityInGround(
            s,
            m_terminationStation);
        Vec3 x_GP = dataPos.terminationPoint_G;

        // Compute lengthening.
        lengthDot += dot(UnitVec3(x_GP - x_GQ), v_GP - v_GQ);
    }

    //--------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------

    // Reference back to the subsystem.
    CableSubsystem* m_subsystem = nullptr;
    CableSpanIndex m_index      = CableSpanIndex::Invalid();

    // TOPOLOGY CACHE (set during realizeTopology())
    CacheEntryIndex m_indexDataInst = CacheEntryIndex::Invalid();
    CacheEntryIndex m_indexDataPos  = CacheEntryIndex::Invalid();
    CacheEntryIndex m_indexDataVel  = CacheEntryIndex::Invalid();

    // Path origin attachment point (a body and a station on the body).
    MobilizedBodyIndex m_originBody = MobilizedBodyIndex::Invalid();
    Vec3 m_originStation{NaN};

    // Path termination attachment point (a body and a station on the body).
    MobilizedBodyIndex m_terminationBody = MobilizedBodyIndex::Invalid();
    Vec3 m_terminationStation{NaN};

    // Path CurveSegments over the obstacles.
    Array_<CurveSegment, ObstacleIndex> m_curveSegments;

    // Path ViaPoints.
    Array_<ViaPoint, ViaPointIndex> m_viaPoints;

    // Path CableSegments separated by via points.
    Array_<CableSegment, CableSegmentIndex> m_cableSegments;

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
    const CableSegment& cableSegment = cable.getCableSegment(
        getCableSegmentIndex());
    const Array_<ObstacleIndex, ObstacleIndex>& obstacleIndexes =
        cableSegment.getObstacleIndexes();

    // If there are no obstacles in this cable segment, there can't be any
    // previous obstacles in contact.
    if (obstacleIndexes.empty()) {
        return ObstacleIndex::Invalid();
    }

    // Find the first obstacle that makes contact before this curve segment
    // in the current cable segment.
    for (ObstacleIndex ix(m_obstacleIndex);
            ix > obstacleIndexes.front(); --ix) {
        ObstacleIndex prevIx(ix - 1);
        const CurveSegment& curveSegment =
            cable.getObstacleCurveSegment(prevIx);
        if (curveSegment.isInContactWithSurface(state)) {
            return prevIx;
        }
    }

    // No active obstacles before this segment.
    return ObstacleIndex::Invalid();
}

ObstacleIndex CurveSegment::findNextObstacleInContactWithCable(
    const State& state) const
{
    const CableSpan::Impl& cable = getCable();
    const CableSegment& cableSegment = cable.getCableSegment(
        getCableSegmentIndex());
    const Array_<ObstacleIndex, ObstacleIndex>& obstacleIndexes =
        cableSegment.getObstacleIndexes();

    // If there are no obstacles in this cable segment, there can't be any
    // next obstacles in contact.
    if (obstacleIndexes.empty()) {
        return ObstacleIndex::Invalid();
    }

    // Find the first active obstacle after this curve segment in the current
    // cable segment.
    for (ObstacleIndex ix(m_obstacleIndex + 1);
            ix <= obstacleIndexes.back(); ++ix) {
        const CurveSegment& curveSegment = cable.getObstacleCurveSegment(ix);
        if (curveSegment.isInContactWithSurface(state)) {
            return ix;
        }
    }

    // No active obstacles after this curve segment.
    return ObstacleIndex::Invalid();
}

Vec3 CurveSegment::findPrevPathPoint_G(const State& state) const
{
    const CableSpan::Impl& cable = getCable();
    const CableSegment& cableSegment = cable.getCableSegment(
        getCableSegmentIndex());
    const Array_<ObstacleIndex, ObstacleIndex>& obstacleIndexes =
        cableSegment.getObstacleIndexes();

    // Check if there is a curve segment preceding this curve segment in the
    // current cable segment.
    if (!obstacleIndexes.empty()) {
        for (ObstacleIndex ix(m_obstacleIndex);
                ix > obstacleIndexes.front(); --ix) {
            ObstacleIndex prevIx(ix - 1);
            const CurveSegment& curveSegment =
                cable.getObstacleCurveSegment(prevIx);
            if (curveSegment.isInContactWithSurface(state)) {
                return curveSegment.calcFinalContactPoint_G(state);
            }
        }
    }

    // Check to see if this cable segment starts with a via point.
    if (cableSegment.getInitialViaPointIndex().isValid()) {
        const ViaPoint& viaPoint = cable.getViaPoint(
            cableSegment.getInitialViaPointIndex());
        return viaPoint.updStation_G(state);
    }

    // There is no via point at the start of this cable segment: the previous
    // point is the path's origin point.
    return getCable()
        .getOriginBody()
        .getBodyTransform(state)
        .shiftFrameStationToBase(cable.getOriginStation());
}

Vec3 CurveSegment::findNextPathPoint_G(const State& state) const
{
    const CableSpan::Impl& cable = getCable();
    const CableSegment& cableSegment = cable.getCableSegment(
        getCableSegmentIndex());
    const Array_<ObstacleIndex, ObstacleIndex>& obstacleIndexes =
        cableSegment.getObstacleIndexes();

    // Check if there is a curve segment following this curve segment in the
    // current cable segment.
    if (!obstacleIndexes.empty()) {
        for (ObstacleIndex ix(m_obstacleIndex + 1);
                ix <= obstacleIndexes.back(); ++ix) {
            const CurveSegment& curveSegment = cable.getObstacleCurveSegment(ix);
            if (curveSegment.isInContactWithSurface(state)) {
                return curveSegment.calcInitialContactPoint_G(state);
            }
        }
    }

    // Check to see if this cable segment ends in a via point.
    if (cableSegment.getFinalViaPointIndex().isValid()) {
        const ViaPoint& viaPoint = cable.getViaPoint(
            cableSegment.getFinalViaPointIndex());
        return viaPoint.updStation_G(state);
    }

    // There is no via point at the end of this cable segment: the next point
    // is the path's termination point.
    return cable.getTerminationBody()
        .getBodyTransform(state)
        .shiftFrameStationToBase(cable.getTerminationStation());
}

//==============================================================================
//                      CABLE FORCE COMPUTATIONS
//==============================================================================
// This section contains all force related computations.
namespace
{

/* Computes the unit force exerted by a CurveSegment on the obstacle body
represented in ground frame. The force can then be obtained by multiplication
with the cable tension. */
void calcUnitForceExertedByCurve(
    const CurveSegment& curve,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CurveSegmentData::Position& dataPos = curve.getDataPos(s);
    const Vec3& x_GB = curve.getMobilizedBody().getBodyOriginLocation(s);

    // Contact point moment arms in ground.
    const Vec3 r_P = dataPos.X_GP.p() - x_GB;
    const Vec3 r_Q = dataPos.X_GQ.p() - x_GB;

    // Tangent directions at contact points in ground.
    const UnitVec3& t_P = getTangent(dataPos.X_GP);
    const UnitVec3& t_Q = getTangent(dataPos.X_GQ);

    // The tangent at the first contact point is along the negative force
    // direction, while at the last contact point the tangent is along the
    // force direction:
    unitForce_G[0] = r_Q % t_Q - r_P % t_P;
    unitForce_G[1] = t_Q - t_P;
}

void calcUnitForceExertedByViaPoint(
    const ViaPoint& viaPoint,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CableSpanData::Position& dataPos = viaPoint.getCable().getDataPos(s);
    ViaPointIndex ix = viaPoint.getIndex();

    const Vec3& x_GB = viaPoint.getMobilizedBody().getBodyOriginLocation(s);
    const UnitVec3& inTangent_G = dataPos.viaPointInTangents_G[ix];
    const UnitVec3& outTangent_G = dataPos.viaPointOutTangents_G[ix];

    // The via point moment arm in ground.
    const Vec3 r_G = viaPoint.getStation_G(s) - x_GB;

    // The incoming direction is along the negative force direction, while the
    // outgoing direction is along the force direction:
    unitForce_G[0] = r_G % outTangent_G - r_G % inTangent_G;
    unitForce_G[1] = outTangent_G - inTangent_G;
}

/* Computes the unit force exerted by the cable at the cable origin, on the
origin body, represented in ground frame. The force can then be obtained by
multiplication with the cable tenson. */
void calcUnitForceAtCableOrigin(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CableSpanData::Position& dataPos = cable.getDataPos(s);

    // Origin contact point moment arm in ground.
    const Vec3 arm_G =
        dataPos.originPoint_G - cable.getOriginBody().getBodyOriginLocation(s);

    // The tangent at the cable origin is along the force direction:
    unitForce_G[0] = arm_G % dataPos.originTangent_G;
    unitForce_G[1] = Vec3(dataPos.originTangent_G);
}

/* Computes the unit force exerted by the cable at the cable termination, on the
termination body, represented in ground frame. The force can then be obtained by
multiplication with the cable tenson. */
void calcUnitForceAtCableTermination(
    const CableSpan::Impl& cable,
    const State& s,
    SpatialVec& unitForce_G)
{
    const CableSpanData::Position& dataPos = cable.getDataPos(s);

    const Vec3 arm_G = dataPos.terminationPoint_G -
                       cable.getTerminationBody().getBodyOriginLocation(s);

    // The tangent at the cable termination is along the negative force
    // direction:
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

    // Forces applied to each via point.
    for (const ViaPoint& viaPoint : m_viaPoints) {
        calcUnitForceExertedByViaPoint(viaPoint, s, unitForce_G);
        viaPoint.getMobilizedBody().applyBodyForce(
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

    for (const ViaPoint& viaPoint : m_viaPoints) {
        calcUnitForceExertedByViaPoint(viaPoint, state, unitForce_G);
        SpatialVec v = viaPoint.getMobilizedBody().getBodyVelocity(state);
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
//                             Path Error Vector and Jacobian
//==============================================================================
/* This section contains helper functions for computing the path error vector
and it's Jacobian to the natural geodesic corrections.

For derivations and background we will refer to the relevant parts in the
original paper: "A fast multi-obstacle muscle wrapping method using natural
geodesic variations" by Andreas Scholz et al (Scholz2015). */
namespace
{

/* Compute the straight line segments in between the obstacles and via points.

Denoted by "e^i" and "e^{i+1} in Scholz2015 Fig 1. */
void calcLineSegments(
    const CableSegment& cableSegment,
    const State& s,
    std::vector<LineSegment>& lines)
{
    lines.clear();
    const CableSpan::Impl& cable = cableSegment.getCable();
    Vec3 prevPathPoint(
        cable.calcInitialCableSegmentPointInGround(s, cableSegment));

    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }
        CurveSegmentData::Position& dataPos = curve.updDataPos(s);

        // Compute the line segment from the previous point to the
        // first curve contact point.
        lines.emplace_back(prevPathPoint, dataPos.X_GP.p());

        // The next line segment will start at the curve's final contact
        // point.
        prevPathPoint = dataPos.X_GQ.p();
    }

    // Compute the last line segment.
    lines.emplace_back(prevPathPoint,
        cable.calcFinalCableSegmentPointInGround(s, cableSegment));
}

/* Helper for computing a single element of the path error vector.
The algorithm for finding the correct path will attempt to drive this path
error to zero.

These are the elements in equation 14 in Scholz2015. */
Real calcPathErrorElement(
    const LineSegment& e,
    const FrenetFrame& X,
    CoordinateAxis axis)
{
    return dot(e.direction, X.R().getAxisUnitVec(axis));
}

/* Compute the path error for each curve segment at the boundary frames, and
stack them in a vector.

If axes = {NormalAxis, BinormalAxis} this gives equation 14 in Scholz2015 */
template <size_t N>
void calcPathErrorVector(
    const CableSegment& cableSegment,
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError)
{
    // Reset path error vector to zero.
    pathError.setToZero();

    // The element in the path error vector to write to.
    int pathErrorIx = -1;
    // Index to a straight line segment.
    int lineIx = 0;

    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        // Path errors are only computed for obstacles that make contact.
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }
        CurveSegmentData::Position& dataPos = curve.updDataPos(s);
        // Compute path error at first contact point.
        for (CoordinateAxis axis : axes) {
            pathError(++pathErrorIx) =
                calcPathErrorElement(lines.at(lineIx), dataPos.X_GP.R(), axis);
        }
        ++lineIx;
        // Compute path error at final contact point.
        for (CoordinateAxis axis : axes) {
            pathError(++pathErrorIx) =
                calcPathErrorElement(lines.at(lineIx), dataPos.X_GQ.R(), axis);
        }
    }
}

} // namespace

//==============================================================================
//                             Path Error Jacobian
//==============================================================================
namespace
{

/* Computes the Jacobian of the path error vector to the
NaturalGeodesicCorrections of all CurveSegments.

When axes = {Normal, Binormal}, this equals equation 53 in Scholz2015. */
template <size_t N>
void calcPathErrorJacobian(
    const CableSegment& cableSegment,
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J)
{
    // Number of free coordinates for a generic geodesic.
    constexpr int NQ = c_GeodesicDOF;

    // Reset the values in the Jacobian.
    J.setToZero();

    // Current indexes to write the elements of the Jacobian to,
    int row = 0;
    int col = 0;

    // Helper for adding a computed block to the Jacobian.
    auto AddBlock = [&](const Vec4& block, int colOffset = 0)
    {
        for (int ix = 0; ix < 4; ++ix) {
            J(row, col + colOffset + ix) += block[ix];
        }
    };

    // Index to a straight line segment.
    int lineIx = 0;
    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }
        CurveSegmentData::Position& dataPos = curve.updDataPos(s);

        // Compute the Jacobian elements of the path error related to the
        // initial contact point.
        for (CoordinateAxis axis : axes) {
            // The path error of which to take the Jacobian is: dot(a_P, e_P)
            const UnitVec3& a_P = dataPos.X_GP.R().getAxisUnitVec(axis);
            const UnitVec3& e_P = lines.at(lineIx).direction;
            const Real l_P      = lines.at(lineIx).length;

            // Compute the variation of the path error.
            const Vec3 dc_dv_P = (a_P - e_P * dot(e_P, a_P)) / l_P;
            const Vec3 dc_dw_P = cross(a_P, e_P);
            // Write to the block diagonal of J the Jacobian of the path error
            // to the CurveSegment's NaturalGeodesicCorrection.
            AddBlock(~dataPos.v_P * dc_dv_P + ~dataPos.w_P * dc_dw_P);

            // Write to the off-diagonal of J the Jacobian of the path
            // error to the previous CurveSegment's NaturalGeodesicCorrection.
            const ObstacleIndex prevObsIx =
                curve.findPrevObstacleInContactWithCable(s);
            if (prevObsIx.isValid()) {
                constexpr int colShift = -NQ; // the off-diagonal block.

                const Vec3 dc_dv_Q = -dc_dv_P;
                const Mat34& v_Q =
                    cable.getObstacleCurveSegment(prevObsIx).updDataPos(s).v_Q;

                AddBlock(~v_Q * dc_dv_Q, colShift);
            }
            // Increment to the next path error element.
            ++row;
        }

        // Compute the Jacobian elements of the path error related to the final
        // contact point.
        ++lineIx; // Increment to the next straight line segment.
        for (CoordinateAxis axis : axes) {
            // The path error of which to take the Jacobian is: dot(a_Q, e_Q)
            const UnitVec3& a = dataPos.X_GQ.R().getAxisUnitVec(axis);
            const UnitVec3& e = lines.at(lineIx).direction;
            const Real l      = lines.at(lineIx).length;

            // Compute the variation of the path error.
            const Vec3 dc_dv_Q = -(a - e * dot(e, a)) / l;
            const Vec3 dc_dw_Q = cross(a, e);
            // Write to the block diagonal of J the Jacobian of the path error
            // to the CurveSegment's NaturalGeodesicCorrection.
            AddBlock(~dataPos.v_Q * dc_dv_Q + ~dataPos.w_Q * dc_dw_Q);

            // Write to the off-diagonal of J the Jacobian of the path
            // error to the next CurveSegment's NaturalGeodesicCorrection.
            const ObstacleIndex nextObsIx =
                curve.findNextObstacleInContactWithCable(s);
            if (nextObsIx.isValid()) {
                constexpr int colShift = NQ; // the off-diagonal block.

                const Vec3 dc_dv_P = -dc_dv_Q;
                const Mat34& v_P =
                    cable.getObstacleCurveSegment(nextObsIx).updDataPos(s).v_P;
                AddBlock(~v_P * dc_dv_P, colShift);
            }
            // Increment to the next path error element.
            ++row;
        }

        col += NQ;
    }
}

/* Computes the gradient of the total cable length to the
NaturalGeodesicCorrections of all CurveSegments.

This function can be derived by viewing the cable as a collection of straight
line segments connecting curve segments. Each curve segment has 4 degrees of
freedom: The NaturalGeodesicCorrections. For each curve segment we can thus
compute the variation of the length of the straight line segment at the initial
Frenet Frame P, the variation of the length of the curve segment, and the
variation of the length of the straight line segment at the final Frenet frame
Q. For one obstacle the gradient is then given as:

g_i = ~e_P * v_P + ~e_Q * v_Q + ~[0, 0, 0, 1]

where P, Q denotes the initial and final frenet frame, e denotes the direction
of the straight line segment, and v denotes the variation of the Frenet frame
position, and g_i is the block of the gradient related to that obstacle.
Stacking of all blocks gives the total gradient. */
void calcLengthGradient(
    const State& state,
    const CableSegment& cableSegment,
    const std::vector<LineSegment>& lines,
    Vector& gradient)
{
    // Reset path error vector to zero.
    gradient.setToZero();

    // The row element in the path error vector to write to.
    int row = -1;
    // Index to a straight line segment.
    int lineIx = 0;

    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        // Only consider obstacles that make contact.
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(state)) {
            continue;
        }
        CurveSegmentData::Position& dataPos = curve.updDataPos(state);

        // Grab the straight line segment's directions.
        const UnitVec3& e_P = lines.at(lineIx).direction;
        const UnitVec3& e_Q = lines.at(lineIx + 1).direction;

        // Compute variation of the total length to the natural geodesic
        // variations of the current CurveSegment.
        Vec4 gradientBlock = ~dataPos.v_P * e_P - ~dataPos.v_Q * e_Q;
        for (int i = 0; i < 4; ++i) {
            gradient(++row) = gradientBlock(i) + (i == 3 ? 1. : 0.);
        }
        ++lineIx;
    }
}

/* Computes the Hessian of the total cable length to the
NaturalGeodesicCorrections of all CurveSegments.

This function can be derived taking the Jacobian of the gradient to the
NaturalGeodesicCorrections of all curve segments. */
void calcLengthHessian(
    const State& state,
    const CableSegment& cableSegment,
    const std::vector<LineSegment>& lines,
    Matrix& hessian)
{
    // Number of free coordinates for a generic geodesic.
    constexpr int NQ = c_GeodesicDOF;

    // Reset the values in the Jacobian.
    Matrix& H = hessian;
    H.setToZero();

    // Current indexes to write the elements of the Jacobian to,
    int row = 0;
    int col = 0;

    // Helper for adding a computed block to the Jacobian.
    auto AddBlock = [&](const Mat44& block, int colOffset = 0)
    {
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                H(row + r, col + colOffset + c) += block(r, c);
            }
        }
    };

    Mat33 I(1.);

    // Helper for constructing a skew symmetric matrix S(v) from a vector v.
    auto skewSymMat = [](const Vec3& v) -> Mat33
    {
        return {
            0., -v(2), v(1),
            v(2), 0., -v(0),
            -v(1), v(0), 0.,
        };
    };

    // Index to a straight line segment.
    int lineIx = 0;
    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(state)) {
            continue;
        }
        CurveSegmentData::Position& dataPos = curve.updDataPos(state);
        const CurveSegmentData::Instance& dataInst = curve.getDataInst(state);

        // This block computes the Hessian of the straight line segment
        // connected to the initial contact point, i.e. at the P-frame.
        {
            const UnitVec3& e = lines.at(lineIx).direction;
            const Real l      = lines.at(lineIx).length;

            const Mat33 E = (I - e * ~e) / l;
            const Mat33 S = skewSymMat(e);

            const Mat34& v_P = dataPos.v_P;
            const Mat34& w_P = dataPos.w_P;

            // The gradient block is given as: g_i = ~e v_i
            // We can use the product rule using the derivatives:
            // de/dq_i = (I - e * ~e) / l * v_i
            // dv_i/dq_j = w_j % v_i
            // Where the latter is true because v_i is a unit vector.
            // The Hessian elements of this block are then given as:
            // H_ij = ~e dv_i/dq_j + ~v_i de/dq_j
            // Which can be written compactly as:
            const Mat44 hblock = ~v_P * E * v_P + ~v_P * S * w_P;
            AddBlock(hblock);

            // Write to the off-diagonal block of H the Jacobian of the
            // gradient to the previous CurveSegment's
            // NaturalGeodesicCorrection.
            const ObstacleIndex prevObsIx =
                curve.findPrevObstacleInContactWithCable(state);
            if (prevObsIx.isValid()) {
                constexpr int colShift = -NQ; // the off-diagonal block.
                const Mat34& v_Q = cable.getObstacleCurveSegment(prevObsIx)
                                       .updDataPos(state)
                                       .v_Q;
                AddBlock(-~v_P * E * v_Q, colShift);
            }
        }

        // This block computes the Hessian of the straight line segment
        // connected to the final contact point, i.e. at the Q-frame.
        ++lineIx; // Increment to the next straight line segment.
        {
            const UnitVec3& e = lines.at(lineIx).direction;
            const Real l      = lines.at(lineIx).length;

            const Mat33 E = (I - e * ~e) / l;
            const Mat33 S = skewSymMat(e);

            const Mat34& v_Q = dataPos.v_Q;
            const Mat34& w_Q = dataPos.w_Q;

            const Real aDot = dataInst.jacobiDot_Q[0];
            const Real rDot = dataInst.jacobiDot_Q[1];

            // If we assume the variation of the jacobi fields is zero, then the
            // columns of v_Q have a constant norm. Then, the Hessian of the
            // straight line segment at the Q-frame is computed analogously as
            // done at the P-frame:
            Mat44 hblock = ~v_Q * E * v_Q - ~v_Q * S * w_Q;

            // In reality the columns of v_Q do not have a constant norm, using
            // the product rule we need to add some additional terms.
            // When computing the variation of the jacobi fields to the natural
            // geodesic variations we need the jacobi covariance: i.e. the
            // second order variation field. We know these terms for the first
            // and last natural geodesic variation, these are simply aDot,
            // rDot. However for the second and third terms these must be
            // computed numerically, similarly to how the jacobi field is
            // computed.
            // Here we assume the effect of these terms to be minimal, and near
            // the optimum these terms will vanish analytically. We therefore
            // estimate the covariance terms as zero. This assumption is
            // verified by the perturbation test done by
            // CableSubsystemTestHelper, which numerically checks that the
            // computed Hessian predicts the variation of the gradient.
            const Real unknown = 0.; // The unknown terms we estimate as zero.
            // Variation of jacobi scalars to each natural geodesic variation:
            const Vec4 da_dq(aDot, unknown, unknown, aDot); // Translational
            const Vec4 dr_dq(rDot, unknown, unknown, rDot); // Rotational
            // All terms related to the covariance are pre-multiplied by this
            // term (e^T * b = dot product of line direction with binormal),
            // which converges to zero near the optimum.
            const Real eTb = dot(e, getBinormal(dataPos.X_GQ));
            // Add the effect of the jacobi field covariance to the hessian
            // block.
            for (int qIx = 0; qIx < c_GeodesicDOF; ++qIx) {
                hblock(1, qIx) += -eTb * da_dq(qIx);
                hblock(2, qIx) += -eTb * dr_dq(qIx);
            }
            AddBlock(hblock);

            // Write to the off-diagonal block of H the Jacobian of the
            // gradient to the next CurveSegment's
            // NaturalGeodesicCorrection.
            const ObstacleIndex nextObsIx =
                curve.findNextObstacleInContactWithCable(state);
            if (nextObsIx.isValid()) {
                constexpr int colShift = NQ; // the off-diagonal block.
                const Mat34& v_P = cable.getObstacleCurveSegment(nextObsIx)
                                       .updDataPos(state)
                                       .v_P;
                AddBlock(-~v_Q * E * v_P, colShift);
            }
        }

        col += NQ;
        row += NQ;
    }
}

} // namespace

//==============================================================================
//                      Solving for NaturalGeodesicCorrections
//==============================================================================

namespace
{

/* Given a correction vector computed from minimizing the cost function,
compute the maximum allowed stepsize along that descending direction. The
allowed step is computed using the orientation variation and the correction
vector. Multiplication of these two gives a rotation vector. The stepsize is
then clamped such that the inf-norm of that rotation vector does not exceed the
maximum allowed angular stepsize. */
void calcCurveCorrectionStepSize(
    const State& s,
    const CurveSegment& curve,
    const NaturalGeodesicCorrection& c,
    Real allowedAngularDisplacement,
    Real& maxAllowedStepSize)
{
    CurveSegmentData::Position& dataPos = curve.updDataPos(s);

    // Helper for computing the maximum angular displacement estimate given the
    // orientation variation at the boundary frame w.
    auto CalcAbsAngularDisplacementEst = [&](const Mat34& w) -> Real
    {
        Real maxAngularDisplacement = 0;
        for (int i = 0; i < w.ncol(); ++i) {
            const Real wAbsMax = max(w.col(i).abs());
            maxAngularDisplacement =
                std::max(maxAngularDisplacement, wAbsMax * c(i));
        }
        return maxAngularDisplacement;
    };

    // Compute the max angular displacement estimate given the angular
    // displacement at the initial and final boundary frame.
    const Real maxAngularDisplacement = std::max(
        CalcAbsAngularDisplacementEst(dataPos.w_P),
        CalcAbsAngularDisplacementEst(dataPos.w_Q));

    // If the estimated angular displacement exceeds the allowed angular
    // displacement, then, limit the stepsize.
    if (maxAngularDisplacement > allowedAngularDisplacement) {
        const Real allowedStepSize =
            std::abs(allowedAngularDisplacement / maxAngularDisplacement);
        maxAllowedStepSize = std::min(maxAllowedStepSize, allowedStepSize);
    }
}

// Helper for computing the stepsize such that a given path correction vector
// does not exceed the max allowed angular displacement at the boundary frames
// of the curve segments.
Real calcPathCorrectionStepSize(
    const State& s,
    const CableSegment& cableSegment,
    const Vector& pathCorrection)
{
    Real stepSize     = 1.;
    int activeCurveIx = 0; // Index of curve in all that are in contact.
    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (!curve.isInContactWithSurface(s)) {
            continue;
        }
        // Each curve segment will evaluate if the stepsize is not too large
        // given the local curvature.
        calcCurveCorrectionStepSize(
            s,
            curve,
            MatrixWorkspace::getCurveCorrection(pathCorrection, activeCurveIx),
            cable.getParameters().solverMaxStepSize,
            stepSize);
        ++activeCurveIx;
    }

    return stepSize;
}

// Helper for computing the length of a cable segment.
// This is the sum of the lengths of all straight line segments, and all curve
// segments comprising the cable segment.
Real calcCableSegmentLength(
    const State& s,
    const CableSegment& cableSegment,
    const std::vector<LineSegment>& lineSegments)
{
    Real cableSegmentLength = 0.;
    // Add length of each line segment.
    for (const LineSegment& line : lineSegments) {
        cableSegmentLength += line.length;
    }
    // Add length of each curve segment.
    const CableSpan::Impl& cable = cableSegment.getCable();
    for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
        const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
        if (curve.isInContactWithSurface(s)) {
            cableSegmentLength += curve.getDataInst(s).arcLength;
        }
    }
    return cableSegmentLength;
};

} // namespace

//------------------------------------------------------------------------------
//                      CableSpan::Impl Cache Computation
//------------------------------------------------------------------------------

void CableSpan::Impl::calcSolverStep(
    const State& s,
    const CableSegment& cableSegment,
    CableSpanAlgorithm algorithm,
    MatrixWorkspace& data) const
{
    // SOLVER STEP 1: Compute the path errors. If the path errors are small the
    // current path is the optimal path, and there is nothing to do. Otherwise
    // we proceed with step 2: computing the corrections to reduce the path
    // errors.

    // The path errors are computed as the misalignment of the straight line
    // segments with the curve segments.
    // Start by computing the straight-line segments of this cable span.
    calcLineSegments(
        cableSegment,
        s,
        data.lineSegments);

    data.maxPathError = 0.; // Reset the max path error field.
    // Only compute the path errors if there are obstacles in contact with the
    // cable.
    if (data.numObstaclesInContact > 0) {
        calcPathErrorVector<1>(
            cableSegment,
            s,
            data.lineSegments,
            {NormalAxis},
            data.normalPathError);
        calcPathErrorVector<1>(
            cableSegment,
            s,
            data.lineSegments,
            {BinormalAxis},
            data.binormalPathError);
        data.maxPathError = std::max(
            data.normalPathError.normInf(),
            data.binormalPathError.normInf());
    }

    // If the path error is small we have converged to the optimal solution,
    // and there is no need to compute the geodesic corrections.
    data.converged = data.maxPathError <= getParameters().smoothnessTolerance;
    if (data.converged) {
        return;
    }

    // SOLVER STEP 2: Compute the path correction step that reduces the path
    // errors.
    SimTK_ASSERT_ALWAYS(
        data.numObstaclesInContact > 0,
        "No obstacles in contact with the cable: unable to compute path corrections");
    // The path corrections can be computed using different algorithms, see
    // CableSpanAlgorithm's documentation for details.
    switch (algorithm) {
    // Compute the path corrections as outlined in Scholz2015.
    case CableSpanAlgorithm::Scholz2015: {
        // Number of path error constraints per curve segment.
        static constexpr int c_NumPathErrorConstraints = 4;
        // Total number of constraints per curve segment.
        static constexpr int c_NumConstraints = c_NumPathErrorConstraints + 1;

        Vector b(data.numObstaclesInContact * c_NumConstraints, 0.);
        Matrix A(
            data.numObstaclesInContact * c_NumConstraints,
            data.numObstaclesInContact * c_GeodesicDOF,
            0.);
        Vector& q = data.pathCorrection;

        // Fill the upper block of A and b with the path error jacobian and
        // vector.
        calcPathErrorVector<2>(
            cableSegment,
            s,
            data.lineSegments,
            {NormalAxis, BinormalAxis},
            b);
        calcPathErrorJacobian<2>(
            cableSegment,
            s,
            data.lineSegments,
            {NormalAxis, BinormalAxis},
            A);

        // The last elements of the b vector contain the change in length
        // of each curve segment, and require it to remain zero. The change in
        // length is the last element of the NaturalGeodesicCorrection vector.
        // The Jacobian is therefore zero everywhere, except for one (=1) at the
        // element of the third NaturalGeodesicCorrection of the corresponding
        // curve segment. Finally we write a weight instead of a one to obtain a
        // weighted least squares. As the weight we take the current maximum
        // path error. This will heavily penalize changing the length when we
        // are far from the optimal solution, and ramp up convergence as we get
        // closer.
        for (int i = 0; i < data.numObstaclesInContact; ++i) {
            // Determine the row and column of the nonzero element in the
            // Jacobian.
            int r = data.numObstaclesInContact * c_NumPathErrorConstraints + i;
            int c = c_GeodesicDOF * (i + 1) - 1;
            // Write the weight that will penalize changing the curve length.
            const Real weight = data.maxPathError;
            A.set(r, c, weight);
        }

        data.factor = A;
        data.factor.solve(b, q);
        q *= -1.;
        break;
    }
    // Compute the path corrections using the Minimal length algorithm.
    case CableSpanAlgorithm::MinimumLength: {
        // This algorithm solves the Lagrangian equations given by:
        // | Q   J^T |   | q |   | -g |
        // | J   0   | * | λ | = | -c |
        // (see CableSpanAlgorithm's documentation).
        const Vector& e = data.normalPathError; // Computed above.
        Matrix& J       = data.normalPathErrorJacobian;
        Vector& lambda  = data.lambda;
        Vector& q       = data.pathCorrection;
        Vector& g       = data.lengthGradient;
        Matrix& H       = data.lengthHessian;
        Matrix& QInv    = data.lengthHessianInverseEstimate;
        // Matrices required for computing the Lagrangian multipliers.
        Matrix& A = data.A;
        Vector& b = data.b;
        // Matrices required for computing QInv from H.
        Vector& sigma = data.singularValues;
        Matrix& U     = data.leftSingularValues;
        Matrix& V     = data.rightSingularValues;

        // Compute the normal path error jacobian.
        calcPathErrorJacobian<1>(cableSegment, s, data.lineSegments,
            {NormalAxis}, J);

        // Compute the gradient and Hessian of the cable length.
        calcLengthGradient(s, cableSegment, data.lineSegments, g);
        calcLengthHessian(s, cableSegment, data.lineSegments, H);

        // Compute QInv:
        // First solve SVD of the Hessian estimate.
        data.svd = 0.5 * (H + ~H);
        data.svd.getSingularValuesAndVectors(sigma, U, V);
        // Now that we have the svd:
        // (H + ~H) / 2 = U * Σ * V
        // we take as the positive definite matrix Q:
        // Q = U * Σ * ~U
        // Note that Q has the same eigen vectors as (H + ~H)/2, as well as
        // the same eigenvalues up to a sign difference.
        //
        // Since we need the inverse of Q, we compute it from the decomposition:
        // Q^-1 = U * Σ^-1 * ~U
        for (int r = 0; r < U.nrow(); ++r) {
            for (int c = 0; c < U.ncol(); ++c) {
                Real elt = 0.;
                for (int i = 0; i < U.ncol(); ++i) {
                    // Add a small weight w to the singular values to make sure
                    // the inverse exists.
                    const Real weight = data.maxPathError;
                    elt += U(r, i) * U(c, i) / (sigma(i) + weight);
                }
                QInv(r, c) = elt;
            }
        }

        // Compute the Lagrange multipliers:
        // (J^T Q^-1 J) * λ = c - J Q^-1 g = c - J d
        //
        // Write as: A * λ = b
        b = e - J * QInv * g;
        A = J * QInv * J.transpose();
        // Solve for the lagrange multipliers.
        data.factor = A;
        data.factor.solve(b, lambda);

        // Finally compute the path corrections.
        q = QInv * (g + J.transpose() * lambda);

        // Use a damped Newton step.
        q *= -0.5;
        break;
    }
    default:
        SimTK_ASSERT(false, "Unknown CableSpanAlgorithm kind");
    }

    // Compute the stepsize along the descending direction.
    const Real stepSize =
        calcPathCorrectionStepSize(s, cableSegment, data.pathCorrection);
    data.pathCorrection *= stepSize;
}

const CableSpanData::Position& CableSpan::Impl::calcDataPos(
    const State& s) const
{
    // Output of this function.
    CableSpanData::Position& dataPos = updDataPos(s);
    dataPos.cableLength = 0.;
    dataPos.smoothness = 0.;
    dataPos.loopIter = 0;
    dataPos.originPoint_G = calcOriginPointInGround(s);
    dataPos.terminationPoint_G = calcTerminationPointInGround(s);

    // This helper function will extract the useful information from the path
    // solver output, and write to the cached CableSpanData.
    auto updDataPosFromSolverResult =
        [&](const State& s,
            const CableSegment& cableSegment,
            const MatrixWorkspace& workspace,
            int solverLoopCount)
    {
        dataPos.cableLength += calcCableSegmentLength(
            s, cableSegment, workspace.lineSegments);
        dataPos.smoothness = std::max(workspace.maxPathError,
            dataPos.smoothness);
        dataPos.loopIter += solverLoopCount;

        // An invalid initial via point index indicates that the segment begins
        // with the cable origin point.
        ViaPointIndex initialViaPointIndex =
            cableSegment.getInitialViaPointIndex();
        if (initialViaPointIndex.isValid()) {
            dataPos.viaPointOutTangents_G[initialViaPointIndex] =
                workspace.lineSegments.front().direction;
        } else {
            dataPos.originTangent_G =
                workspace.lineSegments.front().direction;
        }

        // An invalid final via point index indicates that the segment ends
        // with the cable termination point.
        ViaPointIndex finalViaPointIndex = cableSegment.getFinalViaPointIndex();
        if (finalViaPointIndex.isValid()) {
            dataPos.viaPointInTangents_G[finalViaPointIndex] =
                workspace.lineSegments.back().direction;
        } else {
            dataPos.terminationTangent_G =
                workspace.lineSegments.back().direction;
        }
    };

    // Compute the station positions of all via points. The via point stations
    // do not change during the solver loop, so we can compute them once here.
    for (const ViaPoint& viaPoint : m_viaPoints) {
        viaPoint.calcStation_G(s);
    }

    // Iterate over all cable segments and solve the sub-problem associated
    // with each.
    for (CableSegmentIndex cix(0); cix < getNumCableSegments(); ++cix) {
        const CableSegment& cableSegment = getCableSegment(cix);
        const Array_<ObstacleIndex>& obstacleIndexes =
            cableSegment.getObstacleIndexes();

        // Start the solver loop for the sub-problem associated with this
        // cable segment. This loop calls calcSolverStep each iteration, and
        // stops when converged or the max iterations was reached.
        for (int loopIter = 0;
            loopIter <=
            getParameters().solverMaxIterations; // Correctly handles zero case.
            ++loopIter) {

            // Make sure all curve segments are realized to position stage.
            // This will transform all last computed geodesics to Ground frame,
            // and will update each curve's WrappingStatus.
            for (ObstacleIndex ix : obstacleIndexes) {
                getObstacleCurveSegment(ix).calcDataPos(s);
            }

            // Grab the matrix workspace used by the solver.
            MatrixWorkspace& workspace = updDataInst(s).updOrInsert(
                countActive(s, cableSegment));

            // Compute the path corrections required to reach the optimal path.
            calcSolverStep(s, cableSegment, getParameters().algorithm,
                workspace);

            // Check if we should stop iterating on this sub-problem.
            const bool maxIterationsReached =
            loopIter >= getParameters().solverMaxIterations;
            if (workspace.converged || maxIterationsReached) {
                updDataPosFromSolverResult(s, cableSegment, workspace, loopIter);
                break;
            }

            // Apply the corrections to each CurveSegment to reduce the path
            // error.
            const Vector& pathCorrection = workspace.pathCorrection;
            int activeCurveIx = 0; // Index of curve in all that are in contact.
            for (ObstacleIndex ix : obstacleIndexes) {
                const CurveSegment& curve = getObstacleCurveSegment(ix);
                if (curve.isInContactWithSurface(s)) {
                    curve.applyGeodesicCorrection(
                        s,
                        MatrixWorkspace::getCurveCorrection(
                            pathCorrection,
                            activeCurveIx));
                    ++activeCurveIx;
                }
            }

            // The applied corrections have changed the path: invalidate each
            // curve segment's cache for this cable segment.
            for (ObstacleIndex ix : obstacleIndexes) {
                // Also invalidate non-active segments: They might touchdown
                // again.
                getObstacleCurveSegment(ix).invalidatePosEntry(s);
            }
        }
    }

    return dataPos;
}

//==============================================================================
//                         CableSubsystemTestHelper
//==============================================================================
// Test correctness of the path error Jacobian by applying a small perturbation
// and comparing the effect on the path error. The perturbation is applied in
// the form of a small geodesic correction along each degree of freedom
// subsequently.
// Only curve segments that are in contact with their respective obstacles are
// tested.
void CableSubsystemTestHelper::Impl::runPerturbationTest(
    const State& state,
    const CableSubsystem& subsystem,
    std::ostream& testReport) const
{
    // Result of this perturbation test.
    bool success = true;

    // Helper lambda for comparing two vectors of ObstacleWrappingStatus.
    auto isEqual = [](const std::vector<ObstacleWrappingStatus>& a,
                      const std::vector<ObstacleWrappingStatus>& b) -> bool
    {
        bool out = a.size() == b.size();
        for (size_t i = 0; i < a.size(); ++i) {
            out = out && a.at(i) == b.at(i);
        }
        return out;
    };

    // Helper lambda returning true if any active curve segment has zero length.
    auto anyCurveSegmentHasZeroLength =
        [&](const State& s, const CableSegment& cableSegment) -> bool
    {
        const CableSpan::Impl& cable = cableSegment.getCable();
        for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
            const CurveSegment& curve = cable.getObstacleCurveSegment(ix);
            if (!curve.isInContactWithSurface(s)) {
                continue;
            }
            if (std::abs(curve.getDataInst(s).arcLength) < Eps) {
                return true;
            }
        }
        return false;
    };

    for (CableSpanIndex cableIx(0); cableIx < subsystem.getNumCables();
         ++cableIx) {
        testReport << "START Jacobian Perturbation Test of Cable " << cableIx
                   << "\n";
        const CableSpan::Impl& cable = subsystem.getCable(cableIx).getImpl();

        // Sanity checks on the config parameters:

        // The Jacobian should predict effect on the path error vector up to 1
        // percent, at least.
        SimTK_ERRCHK1_ALWAYS(
            m_perturbationTestParameters.tolerance <= 1.,
            "CableSubsystemTestHelper::Impl::runPerturbationTest",
            "Perturbation test tolerance (=%e) should be less or equal to 1",
            m_perturbationTestParameters.tolerance);
        // Assertions below are ball park estimates, but what we want to avoid
        // is testing a very small perturbation on a very inaccurate curve
        // segment. That won't work.
        SimTK_ERRCHK2_ALWAYS(
            cable.getParameters().geodesicIntegratorAccuracy <
                m_perturbationTestParameters.perturbation * 1e-3,
            "CableSubsystemTestHelper::Impl::runPerturbationTest",
            "Curve segment accuracy (=%e) should be smaller than the applied perturbation (=%e) times 1e-3",
            cable.getParameters().geodesicIntegratorAccuracy,
            m_perturbationTestParameters.perturbation);
        SimTK_ERRCHK2_ALWAYS(
            m_perturbationTestParameters.perturbation <=
                m_perturbationTestParameters.tolerance * 1e-4,
            "CableSubsystemTestHelper::Impl::runPerturbationTest",
            "Applied perturbation (=%e) should be less or equal to the assertion tolerance (=%e) times 1e-4.",
            m_perturbationTestParameters.perturbation,
            m_perturbationTestParameters.tolerance);

        // Run the perturbation test for each cable segment.
        for (CableSegmentIndex ix(0); ix < cable.getNumCableSegments(); ++ix) {
            testReport << "Testing CableSegment " << ix + 1 << " of "
                       << cable.getNumCableSegments() << "\n";

            // Get the cable segment.
            const CableSegment& cableSegment = cable.getCableSegment(ix);
            const int nActive = cable.countActive(state, cableSegment);

            if (nActive == 0) {
                testReport << "CableSegment has no active segments:";
                testReport << "Skipping perturbation test\n";
                continue;
            }

            if (anyCurveSegmentHasZeroLength(state, cableSegment)) {
                testReport << "CableSegment has zero-length curve segments:";
                testReport << "Skipping perturbation test\n";
                continue;
            }

            // Apply a small perturbation to each of the curve segments in the
            // form of a geodesic correction. There are 4 DOF for each geodesic,
            // and just to be sure we apply a negative and positive perturbation
            // along each DOF of each geodesic.
            Real perturbation = m_perturbationTestParameters.perturbation;
            for (int i = 0; i < c_GeodesicDOF * 2 * nActive + 1; ++i) {
                // We do not want to mess with the actual state, so we make a
                // copy.
                const State sCopy = state;
                cableSegment.getSubsystem().getMultibodySystem().realize(
                    sCopy, Stage::Position);

                // Trigger realizing position level cache, resetting the
                // configuration.
                const CableSpanData::Position& dataPos =
                    cableSegment.getCable().getDataPos(sCopy);

                // The wrapping status of each obstacle is reset after copying
                // the state. The matrix sizes are no longer compatible if this
                // results in a change in wrapping status (if we are touching
                // down during this realization for example). We can only test
                // the Jacobian if the wrapping status of all obstacles remains
                // the same.
                if (cable.countActive(sCopy, cableSegment) != nActive) {
                    testReport <<
                        "Cable wrapping status changed after cloning the state: "
                        "Skipping perturbation test\n";
                    break;
                }

                // Compute the path error vector and Jacobian, before the
                // applied perturbation.
                std::vector<LineSegment> lineSegments(nActive + 1);

                Vector pathError(nActive * c_GeodesicDOF, NaN);
                Matrix pathErrorJacobian(
                    nActive * c_GeodesicDOF,
                    nActive * c_GeodesicDOF,
                    NaN);

                Real length = NaN;
                Vector lengthGradient(nActive * c_GeodesicDOF, NaN);
                Matrix lengthHessian(
                    nActive * c_GeodesicDOF,
                    nActive * c_GeodesicDOF,
                    NaN);

                calcLineSegments(
                    cableSegment,
                    sCopy,
                    lineSegments);

                calcPathErrorVector<2>(
                    cableSegment,
                    sCopy,
                    lineSegments,
                    {NormalAxis, BinormalAxis},
                    pathError);
                calcPathErrorJacobian<2>(
                    cableSegment,
                    sCopy,
                    lineSegments,
                    {NormalAxis, BinormalAxis},
                    pathErrorJacobian);

                length = calcCableSegmentLength(sCopy, cableSegment,
                    lineSegments);
                calcLengthGradient(sCopy, cableSegment, lineSegments,
                    lengthGradient);
                calcLengthHessian(sCopy, cableSegment, lineSegments,
                    lengthHessian);

                // Define the perturbation we will use for testing the Jacobian.
                Vector pathCorrection(nActive * c_GeodesicDOF, 0.);
                if (i != 0) {
                    perturbation *= -1.;
                    pathCorrection.set((i - 1) / 2, perturbation);
                }

                // Use the Jacobian to predict the perturbed path error vector.
                const Vector predictedPathError =
                    pathErrorJacobian * pathCorrection + pathError;
                const Real predictedLength =
                    ~lengthGradient * pathCorrection + length;
                const Vector predictedLengthGradient =
                    lengthHessian * pathCorrection + lengthGradient;

                // Store the wrapping status before applying the perturbation.
                std::vector<ObstacleWrappingStatus> prevWrappingStatus;
                for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                    const CurveSegment& curve =
                        cable.getObstacleCurveSegment(ix);
                    prevWrappingStatus.push_back(
                        curve.getDataInst(sCopy).wrappingStatus);
                }

                // Apply the perturbation.
                int activeCurveIx = 0; // Index of curves that are in contact.
                for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                    const CurveSegment& curve =
                        cable.getObstacleCurveSegment(ix);
                    if (curve.isInContactWithSurface(sCopy)) {
                        curve.applyGeodesicCorrection(
                            sCopy,
                            MatrixWorkspace::getCurveCorrection(
                                pathCorrection,
                                activeCurveIx));
                        ++activeCurveIx;
                    }
                }

                for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                    const CurveSegment& curve =
                        cable.getObstacleCurveSegment(ix);
                    curve.invalidatePosEntry(sCopy);
                }

                for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                    const CurveSegment& curve =
                        cable.getObstacleCurveSegment(ix);
                    curve.calcDataPos(sCopy);
                }

                // Compute the path error vector after having applied the
                // perturbation.
                calcLineSegments(
                    cableSegment,
                    sCopy,
                    lineSegments);

                calcPathErrorVector<2>(
                    cableSegment,
                    sCopy,
                    lineSegments,
                    {NormalAxis, BinormalAxis},
                    pathError);
                calcPathErrorJacobian<2>(
                    cableSegment,
                    sCopy,
                    lineSegments,
                    {NormalAxis, BinormalAxis},
                    pathErrorJacobian);

                length = calcCableSegmentLength(sCopy, cableSegment,
                    lineSegments);
                calcLengthGradient(sCopy, cableSegment, lineSegments,
                    lengthGradient);
                calcLengthHessian(sCopy, cableSegment, lineSegments,
                    lengthHessian);

                // The curve lengths are clipped at zero, which distorts the
                // perturbation test. We simply skip these edge cases.
                if (anyCurveSegmentHasZeroLength(sCopy, cableSegment)) {
                    testReport << "CableSegment has zero-length curve segments after perturbation:";
                    testReport << "Skipping perturbation test\n";
                    continue;
                }

                // We can only use the path error Jacobian if the small
                // perturbation did not trigger any liftoff or touchdown on any
                // obstacles. If any CurveSegment's WrappingStatus has changed,
                // we will not continue with the test. Since the perturbation is
                // small, this is unlikely to happen often.
                std::vector<ObstacleWrappingStatus> nextWrappingStatus;
                for (ObstacleIndex ix : cableSegment.getObstacleIndexes()) {
                    const CurveSegment& curve =
                        cable.getObstacleCurveSegment(ix);
                    nextWrappingStatus.push_back(
                        curve.getDataInst(sCopy).wrappingStatus);
                }
                if (!isEqual(prevWrappingStatus, nextWrappingStatus)) {
                    testReport << "Wrapping status changed: Stopping test\n";
                    break;
                }

                // Tolerance was given as a percentage, multiply by 0.01 to get the
                // correct scale.
                const Real tolerance =
                    m_perturbationTestParameters.tolerance * 1e-2;
                // Evaluate the path error jacobian.
                {
                    const Real predictionError =
                        (predictedPathError - pathError).normInf();
                    bool passedTest =
                        std::abs(predictionError / perturbation) <= tolerance;
                    if (!passedTest) {
                        testReport << "FAILED path error Jacobian perturbation test for correction = "
                                   << pathCorrection << "\n";
                        testReport << "path error after perturbation:\n";
                        testReport << "    Got        : " << pathError << "\n";
                        testReport << "    Predicted  : " << predictedPathError
                                   << "\n";
                        testReport << "    Difference : "
                                   << predictionError / perturbation
                                   << "\n";
                        testReport << "    Max diff   : "
                                   << predictionError / perturbation
                                   << " <= " << tolerance << " = eps\n";
                    } else {
                        testReport << "PASSED perturbation test for correction = "
                                   << pathCorrection << "\n";
                        testReport << " ( max diff = "
                                   << predictionError / perturbation
                                   << " <= " << tolerance << " = eps )\n";
                    }
                    success = success && passedTest;
                }
                // Evaluate the length gradient.
                {
                    const Real predictionError =
                        std::abs(predictedLength - length);
                    bool passedTest =
                        std::abs(predictionError / perturbation) <= tolerance;
                    if (!passedTest) {
                        testReport << "FAILED length gradient perturbation test for correction = "
                                    << pathCorrection << "\n";
                        testReport << "length after perturbation:\n";
                        testReport << "    Got        : " << length << "\n";
                        testReport << "    Predicted  : " << predictedLength
                                   << "\n";
                        testReport << "    Difference : "
                                   << predictionError / perturbation
                                   << "\n";
                        testReport << "    Tolerance  : " << tolerance << "\n";
                    } else {
                        testReport << "PASSED length perturbation test for correction = "
                                   << pathCorrection << "\n";
                        testReport << " ( max diff = "
                                   << std::abs(predictionError / perturbation)
                                   << " <= " << tolerance << " = eps )\n";
                    }
                    success = success && passedTest;
                }
                // Evaluate the length Hessian.
                {
                    const Real predictionError =
                        (predictedLengthGradient - lengthGradient).normInf();
                    bool passedTest =
                        std::abs(predictionError / perturbation) <= tolerance;
                    if (!passedTest) {
                        testReport << "FAILED length Hessian perturbation test for correction = "
                                   << pathCorrection << "\n";
                        testReport << "length gradient after perturbation:\n";
                        testReport << "    Got        : "
                                   << lengthGradient << "\n";
                        testReport << "    Predicted  : "
                                   << predictedLengthGradient << "\n";
                        testReport << "    Difference : "
                                   << predictionError / perturbation
                                   << "\n";
                        testReport << "    Max diff   : "
                                   << std::abs(predictionError / perturbation)
                                   << " <= " << tolerance << " = eps\n";
                    } else {
                        testReport << "PASSED length Hessian perturbation test for correction = "
                                   << pathCorrection << "\n";
                        testReport << " ( max diff = "
                                   << std::abs(predictionError / perturbation)
                                   << " <= " << tolerance << " = eps )\n";
                    }
                    success = success && passedTest;
                }
            }
        }
    }

    // Flush all info to the report, in case we throw an exception.
    if (success) {
        testReport << "PASSED TEST: Jacobian perturbation test" << std::endl;
    } else {
        testReport << "FAILED TEST: Jacobian perturbation test" << std::endl;
    }
    SimTK_ASSERT_ALWAYS(success, "Perturbation test failed");
}

//==============================================================================
//                      Resampling geodesic path points
//==============================================================================
/* This section contains functions for resampling a computed geodesic. */

namespace
{

/* Bound for the distance between two knots when interpolation will switch
from a Hermite interpolation to a linear interpolation to retain stability. */
constexpr Real s_SMALL_KNOT_DISTANCE_BOUND = 1e-6;

// Perform Hermite interpolation given the knot = [(y0, y0Dot) @ x0] and
// knot = [(y1, y1Dot) @ x1], to compute y at given x. If the distance dx=x1-x0
// becomes too small this function falls back to a linear interpolation.
Vec3 calcHermiteInterpolation(
    Real x0,
    const Vec3& y0,
    const Vec3& y0Dot,
    Real x1,
    const Vec3& y1,
    const Vec3& y1Dot,
    Real x)
{
    SimTK_ASSERT2(
        x0 <= x,
        "Hermite interpolation domain error: "
        "x0 = %e must be smaller or equal than x = %e",
        x0,
        x);
    SimTK_ASSERT2(
        x <= x1,
        "Hermite interpolation domain error: "
        "x = %e must be smaller or equal than x1 = %e",
        x,
        x1);

    // Distance between the knots:
    const Real dx = x1 - x0;
    // For very very small distance between the knots just return the average of
    // the values.
    if (dx < Eps) {
        return (y0 + y1) / 2.;
    }

    // Normalize the point of interpolation:
    const Real s = (x - x0) / dx; // Range of s is [0, 1].
    // For small distance between the knots return the linear interpolation
    // between the values.
    if (dx < s_SMALL_KNOT_DISTANCE_BOUND) {
        return y0 + (y1 - y0) * s;
    }

    // Otherwise do Hermite interpolation.
    // Helper constants:
    const Vec3 dy = y1 - y0;
    const Vec3 m  = 0.5 * dx * (y0Dot + y1Dot);
    // Compute coefficients of polynomial:
    // y = c0 + c1 * s + c2 * s^2 + c3 * s^3
    const Vec3& c0 = y0;
    const Vec3& c1 = y0Dot * dx;
    const Vec3 c3  = 2. * (m - dy);
    const Vec3 c2  = 3. * dy - 2. * m - c1;
    // Evaluate polynomial at s.
    return c0 + s * (c1 + s * (c2 + s * c3));
}

// Resample the geodesic at equal length intervals, and write the computed
// points to the provided sink.
void calcResampledGeodesicPoints(
    const std::vector<ContactGeometry::GeodesicKnotPoint>& geodesic,
    int numSamples,
    const std::function<void(const Vec3& interpolatedPoint_S)>& sink)
{
    // Some sanity checks.
    SimTK_ASSERT(
        !geodesic.empty(),
        "Resampling of geodesic failed: Provided geodesic is empty.");
    SimTK_ASSERT(
        geodesic.front().arcLength == 0.,
        "Resampling of geodesic failed: First frame must be at arcLength = zero");
    SimTK_ASSERT(
        geodesic.back().arcLength >= 0.,
        "Resampling of geodesic failed: Last frame must be at arcLength >= zero");
    SimTK_ASSERT1(
        numSamples >= 2,
        "Resampling of geodesic failed: Minimum number of samples for resampling is 2 (got %i)",
        numSamples);

    // If there is but one sample (arcLength = 0), write that sample.
    if (geodesic.size() == 1) {
        for (int i = 0; i < numSamples; ++i) {
            sink(geodesic.front().point);
        }
        return;
    }

    // The resampled geodesic will capture the first and last value exactly.
    // Write the first value:
    sink(geodesic.front().point);

    // Since the first and last point do not need interpolation, compute the
    // remaining n=(numSamples-2) points using interpolation:
    auto it                   = geodesic.begin();
    const Real finalArcLength = geodesic.back().arcLength;
    for (int i = 1; i < numSamples - 1; ++i) {
        // Length at the current interpolation point.
        const Real lenghtAtInterpolationPoint =
            finalArcLength *
            (static_cast<Real>(i) / static_cast<Real>(numSamples - 1));

        // Find the two knots (a, b) of the geodesic such that the
        // arcLength of the interpolation point l lies between them.
        // i.e. find: a.arcLength <= l <= b.arcLength
        auto b = std::lower_bound(
            it,
            geodesic.end(),
            lenghtAtInterpolationPoint,
            [&](const ContactGeometry::GeodesicKnotPoint& y, Real x)
            { return y.arcLength < x; });
        auto a = b == geodesic.begin() ? b : b - 1;

        SimTK_ASSERT(
            b != geodesic.end(),
            "Attempting to access out-of-range element");

        // Interpolate between knots (a, b).
        const Vec3 interpolatedPoint = calcHermiteInterpolation(
            a->arcLength,
            a->point,
            a->tangent,
            b->arcLength,
            b->point,
            b->tangent,
            lenghtAtInterpolationPoint);
        // Write interpolated point to the output.
        sink(interpolatedPoint);

        // Update the search range.
        it = a;
    }

    // Write the last value:
    sink(geodesic.back().point);
}

} // namespace

void CurveSegment::calcResampledPoints(
    const State& state,
    int numSamples,
    const std::function<void(Vec3 point_G)>& sink) const
{
    const CurveSegmentData::Instance& dataInst = getDataInst(state);
    const CurveSegmentData::Position& dataPos  = getDataPos(state);
    const ContactGeometry& geometry            = getContactGeometry();

    // If the obstacle's surface is analytic we compute the geodesic with the
    // desired granularity.
    if (geometry.isAnalyticFormAvailable()) {
        geometry.shootGeodesicInDirectionAnalytically(
            dataInst.X_SP.p(),
            getTangent(dataInst.X_SP),
            dataInst.arcLength,
            numSamples,
            // Transform points to ground frame and write to the output sink.
            [&](const ContactGeometry::GeodesicKnotPoint& q)
            {
                const Vec3 point_G =
                    dataPos.X_GS.shiftFrameStationToBase(q.point);
                sink(point_G);
            });
        return;
    }

    // If the obstacle's surface is not analytic we will resample the points
    // from the integrator at equal length intervals using Hermite
    // interpolation.
    calcResampledGeodesicPoints(
        dataInst.geodesicKnotPoints,
        numSamples,
        // Transform points to ground frame and write to the output sink.
        [&](const Vec3& x_S)
        {
            const Vec3 x_G = dataPos.X_GS.shiftFrameStationToBase(x_S);
            sink(x_G);
        });
}

void CableSpan::Impl::calcCurveSegmentResampledPoints(
    const State& state,
    CableSpanObstacleIndex ix,
    int numSamples,
    const std::function<void(Vec3 point_G)>& sink) const
{
    realizePosition(state);
    const CurveSegment& curve = getObstacleCurveSegment(ix);
    SimTK_ERRCHK1_ALWAYS(
        curve.isInContactWithSurface(state),
        "CableSpan::calcCurveSegmentResampledPoints",
        "Obstacle %i is not in contact with cable",
        ix);
    SimTK_ERRCHK1_ALWAYS(
        numSamples >= 2,
        "CableSpan::calcCurveSegmentResampledPoints",
        "Number of resampling points must be larger or equal to 2, but got %i",
        numSamples);
    curve.calcResampledPoints(state, numSamples, sink);
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

    CacheEntryIndex indexDataInst = allocateCacheEntry(
        state,
        Stage::Instance,
        Stage::Infinity,
        new Value<CableSpanData::Instance>());

    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        CableSpan& cable = mutableThis->updCable(ix);
        cable.updImpl().realizeTopology(state, indexDataInst);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemModelImpl(State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeModel(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemInstanceImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeInstance(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemTimeImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeTime(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemPositionImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizePosition(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemVelocityImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeVelocity(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemDynamicsImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeDynamics(state);
    }

    return 0;
}

int
CableSubsystem::Impl::realizeSubsystemAccelerationImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeAcceleration(state);
    }

    return 0;
}

int CableSubsystem::Impl::realizeSubsystemReportImpl(const State& state) const
{
    for (CableSpanIndex ix(0); ix < cables.size(); ++ix) {
        getCable(ix).getImpl().realizeReport(state);
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
    mbs.setCableSubsystem(*this);
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

CableSpan::CableSpan() : m_impl(new Impl())
{}

CableSpan::CableSpan(
    CableSubsystem& subsystem,
    MobilizedBodyIndex originBody,
    const Vec3& originStation,
    MobilizedBodyIndex terminationBody,
    const Vec3& terminationStation) :
    m_impl(new Impl(
        originBody,
        originStation,
        terminationBody,
        terminationStation))
{
    CableSpanIndex ix = subsystem.updImpl().adoptCable(*this);
    updImpl().setSubsystem(subsystem, ix);
    updImpl().createFirstCableSegment();
}

MobilizedBodyIndex CableSpan::getOriginBodyIndex() const
{
    return m_impl->getOriginBodyIndex();
}

void CableSpan::setOriginBodyIndex(MobilizedBodyIndex originBody)
{
    m_impl->setOriginBodyIndex(originBody);
}

MobilizedBodyIndex CableSpan::getTerminationBodyIndex() const
{
    return m_impl->getTerminationBodyIndex();
}

void CableSpan::setTerminationBodyIndex(MobilizedBodyIndex terminationBody)
{
    m_impl->setTerminationBodyIndex(terminationBody);
}

Vec3 CableSpan::getOriginStation() const
{
    return m_impl->getOriginStation();
}

void CableSpan::setOriginStation(const Vec3& originStation)
{
    m_impl->setOriginStation(originStation);
}

Vec3 CableSpan::getTerminationStation() const
{
    return m_impl->getTerminationStation();
}

void CableSpan::setTerminationStation(const Vec3& terminationStation)
{
    m_impl->setTerminationStation(terminationStation);
}

void CableSpan::calcOriginUnitForce(
    const State& state,
    SpatialVec& unitForce_G) const
{
    getImpl().realizePosition(state);
    calcUnitForceAtCableOrigin(getImpl(), state, unitForce_G);
}

void CableSpan::calcTerminationUnitForce(
    const State& state,
    SpatialVec& unitForce_G) const
{
    getImpl().realizePosition(state);
    calcUnitForceAtCableTermination(getImpl(), state, unitForce_G);
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

ViaPointIndex CableSpan::addViaPoint(
    MobilizedBodyIndex viaPointBody,
    const Vec3& station_B)
{
    return updImpl().addViaPoint(viaPointBody, station_B);
}

int CableSpan::getNumViaPoints() const
{
    return getImpl().getNumViaPoints();
}

CableSpanIndex CableSpan::getIndex() const
{
    return getImpl().getIndex();
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

const Transform& CableSpan::getObstacleTransformSurfaceToBody(
    ObstacleIndex ix) const
{
    return getImpl().getObstacleCurveSegment(ix).getTransformSurfaceToBody();
}
void CableSpan::setObstacleTransformSurfaceToBody(
    ObstacleIndex ix,
    const Transform& X_BS)
{
    updImpl().updObstacleCurveSegment(ix).setTransformSurfaceToBody(X_BS);
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

const MobilizedBodyIndex& CableSpan::getViaPointMobilizedBodyIndex(
    ViaPointIndex ix) const
{
    return getImpl().getViaPoint(ix).getMobilizedBodyIndex();
}

void CableSpan::setViaPointMobilizedBodyIndex(
    ViaPointIndex ix,
    MobilizedBodyIndex body)
{
    updImpl().updViaPoint(ix).setMobilizedBodyIndex(body);
}

Vec3 CableSpan::getViaPointStation(ViaPointIndex ix) const
{
    return getImpl().getViaPoint(ix).getStation_B();
}

void CableSpan::setViaPointStation(
    ViaPointIndex ix,
    const Vec3& station_B)
{
    updImpl().updViaPoint(ix).setStation_B(station_B);
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

void CableSpan::setAlgorithm(CableSpanAlgorithm algorithm)
{
    updImpl().updParameters().algorithm = algorithm;
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

void CableSpan::calcCurveSegmentResampledPoints(
    const State& state,
    CableSpanObstacleIndex ix,
    int numSamples,
    const std::function<void(Vec3 point_G)>& sink) const
{
    getImpl().calcCurveSegmentResampledPoints(state, ix, numSamples, sink);
}

void CableSpan::calcCurveSegmentUnitForce(
    const State& state,
    CableSpanObstacleIndex ix,
    SpatialVec& unitForce_G) const
{
    getImpl().realizePosition(state);
    const auto& curve = getImpl().getObstacleCurveSegment(ix);
    calcUnitForceExertedByCurve(curve, state, unitForce_G);
}

Vec3 CableSpan::calcViaPointLocation(
    const State& state,
    CableSpanViaPointIndex ix) const
{
    return getImpl().getViaPoint(ix).getStation_G(state);
}

void CableSpan::calcViaPointUnitForce(
    const State& state,
    CableSpanViaPointIndex ix,
    SpatialVec& unitForce_G) const
{
    getImpl().realizePosition(state);
    const auto& viaPoint = getImpl().getViaPoint(ix);
    calcUnitForceExertedByViaPoint(viaPoint, state, unitForce_G);
}
