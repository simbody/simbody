#ifndef SimTK_SIMBODY_CABLE_SPAN_H_
#define SimTK_SIMBODY_CABLE_SPAN_H_

/* -------------------------------------------------------------------------- *
 *                           Simbody(tm)                                      *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Pepijn van den Bos                                                *
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

#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/common.h"

namespace SimTK
{

SimTK_DEFINE_UNIQUE_INDEX_TYPE(CurveSegmentIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanIndex);

class MultibodySystem;
class CableSubsystem;
class CableSubsystemTestHelper;
class CableSpan;

//==============================================================================
// A curved cable segment that is part of a `CableSpan`, and wraps over an
// obstacle surface.
class SimTK_SIMBODY_EXPORT CurveSegment final
{
private:
    CurveSegment() = default;

public:
    CurveSegment(const CurveSegment&)                = default;
    CurveSegment& operator=(const CurveSegment&)     = default;
    CurveSegment(CurveSegment&&) noexcept            = default;
    CurveSegment& operator=(CurveSegment&&) noexcept = default;
    ~CurveSegment()                                  = default;

    /** Construct a CurveSegment representing a segment of a CableSpan that
    wraps over an obstacle surface (i.e. a ContactGeometry).
    @param cable The cable this segment belongs to.
    @param mobod The body that the contact geometry is rigidly attached to.
    @param X_BS Transform specifying the location and orientation of the
    contact geometry's origin frame with respect to the mobilized body.
    @param geometry The contact geometry over which this segment wraps.
    @param initialContactPointHint A guess of the contact point of the cable
    span and the contact geometry to compute the initial cable path. This point
    is defined in the local contact geometry's frame. The point will be used as
    a starting point when computing the initial cable path. As such, it does
    not have to lie on the contact geometry's surface, nor does it have to
    belong to a valid cable path.*/
    CurveSegment(
        CableSpan cable,
        MobilizedBodyIndex body,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 initialContactPointHint);

    /** Wrapping status of this segment in relation to the contact geometry. */
    enum class WrappingStatus
    {
        InContactWithSurface,
        LiftedFromSurface,
        Disabled,
        // TODO WrappingStatus::InitialGuess is used to initialize the path. Instead this should be done using the first Stage realizations.
        InitialGuess,
    };

    /** Get the ContactGeometry that this segment wraps over. */
    const ContactGeometry& getContactGeometry() const;

    /** Get the MobilizedBody that the contact geometry is rigidly attached to.
     */
    const Mobod& getMobilizedBody() const;

    /** Get the transform representing the orientation and postion of the
    contact geometry's origin with respect to the body fixed frame. */
    const Transform& getXformSurfaceToBody() const;

    /** Set the transform representing the orientation and postion of the
    contact geometry's origin with respect to the body fixed frame. */
    void setXformSurfaceToBody(Transform X_BS);

    /** Get the length of this segment.
    The system must be realiezd to Stage::Position. */
    Real getSegmentLength(const State& state) const;

    /** Get the frenet frame at the start (first contact point) of this curve
    segment.
    TODO describe the frame axes.
    The system must be realiezd to Stage::Position. */
    const Transform& getFrenetFrameStart(const State& state) const;

    /** Get the frenet frame at the end (last contact point) of this curve
    segment.
    TODO describe the frame axes.
    The system must be realiezd to Stage::Position. */
    const Transform& getFrenetFrameEnd(const State& state) const;

    /** Get the wrapping status of this segment.
    The system must be realiezd to Stage::Position. */
    WrappingStatus getStatus(const State& state) const;

    /** Get the number of steps taken by the GeodesicIntegrator to compute this
    segment during the last realization. The system must be realiezd to
    Stage::Position. */
    int getNumberOfIntegratorStepsTaken(const State& state);

    /** Get the initial step size that the GeodesicIntegrator will use for the
    next path computation. The system must be realiezd to Stage::Position. */
    Real getInitialIntegratorStepSize(const State& state);

    // TODO useful?
    /* void setDisabled(const State& state) const; */
    /* void setEnabled(const State& state) const; */

    /** Comute resampled points at equal length intervals along the curve, in
    ground frame.
    The system must be realiezd to Stage::Position.

    @param state State of the system.
    @param points_G The output buffer to which the points are written.
    @param Controls the number of points that are computed. These points are
    resampled from the computed curve at equal length intervals. If the curve
    length is zero, this parameter is ignored, and a single point is written.
    If the curve is not active (Status::Disabled or Status::Lifted), no points
    are written. Note that if nPoints=1, and the curve length is not zero, an
    exception is thrown.
    @return The number of points written. */
    int calcPathPoints(
        const State& state,
        std::function<void(Vec3 point_G)>& sink,
        int nPoints) const;
    int calcPathPointsAndTangents(
        const State& state,
        std::function<void(Vec3 point_G, UnitVec3 tangent_G)>& sink,
        int nSamples = 0) const;

//------------------------------------------------------------------------------
    class Impl;

private:
    friend CableSpan;
    const Impl& getImpl() const
    {
        return *m_Impl;
    }
    Impl& updImpl()
    {
        return *m_Impl;
    }

    std::shared_ptr<Impl> m_Impl = nullptr;
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                          CABLE SPAN
//==============================================================================
// Representation of a cable spanning between two poins.
// The cable can adopt surface obstacles that must be wrapped over.
//
// The CableSpan can be seen as consisting of alternating LineSegments and
// CurveSegments, starting and ending with a straight LineSegment, and a
// CurveSegment for each obstacle it comes into contact with.
//
// NOTE: The interaction with the obstacles is **ordered**. The cable can wrap
// over the obstacles in the order that they are added. This means that the
// cable can not wrap over the first obstacle twice, even though it might
// spatially intersect it twice.
class SimTK_SIMBODY_EXPORT CableSpan final
{
public:

    // Representation of the cable segments that do not lie on a surface.
    struct LineSegment final
    {
        LineSegment() = default;

        LineSegment(Vec3 a, Vec3 b) : length((b - a).norm()), direction((b - a) / length)
        {}

        // TODO rename to length.
        Real length = NaN;
        // TODO rename to direction.
        UnitVec3 direction{NaN, NaN, NaN};
    };


    CableSpan()                                = default;
    ~CableSpan()                               = default;
    CableSpan(const CableSpan&)                = default;
    CableSpan& operator=(const CableSpan&)     = default;
    CableSpan(CableSpan&&) noexcept            = default;
    CableSpan& operator=(CableSpan&&) noexcept = default;

//------------------------------------------------------------------------------

    CableSpan(
        CableSubsystem& subsystem,
        MobilizedBodyIndex originBody,
        const Vec3& defaultOriginPoint,
        MobilizedBodyIndex terminationBody,
        const Vec3& defaultTerminationPoint);

//------------------------------------------------------------------------------

    void addSurfaceObstacle(
        MobilizedBodyIndex mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint = {1., 0., 0.});

    int getNumSurfaceObstacles() const;

//------------------------------------------------------------------------------

    int getNumCurveSegments() const;

    const CurveSegment& getCurveSegment(CurveSegmentIndex ix) const;

    // TODO rename to constraint tolerance
    Real getSurfaceConstraintTolerance() const;
    void setSurfaceConstraintTolerance(Real tolerance);

    // Maximum number of solver iterations for projecting the geodesic state to
    // the surface during shooting of a curve segment over the obstacle's
    // surface.
    // TODO this is not connected to anything currently...
    int getSurfaceProjectionMaxIter() const;
    void setSurfaceProjectionMaxIter(int maxIter);

    // TODO merge with surface constraint tolerance above?
    Real getIntegratorAccuracy() const;
    void setIntegratorAccuracy(Real accuracy);

    // Maximum number of solver iterations for finding the optimal path.
    int getSolverMaxIter() const;
    void setSolverMaxIter(int maxIter);

    // TODO merge with surface constraint tolerance above?
    Real getPathErrorAccuracy() const;
    void setPathErrorAccuracy(Real accuracy);

    // TODO advanced parameter: Max allowed step during Newton iteration in angles. Using the local curvatures the max angular step is converted to a max allowed translation.
    Real getMaxRadialStepInDegrees() const;
    void setMaxRadialStepInDegrees(Real maxStepInDeg);


    Real getLength(const State& state) const;

    Real getLengthDot(const State& state) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    Real calcCablePower(const State& state, Real tension) const;

    // Calls `CurveSegment::calcPoints` on each active curve segment, see
    // description there.
    int calcPathPoints(
        const State& state,
        std::function<void(Vec3 point_G)>& sink,
        int nPointsPerCurveSegment = 0) const;


    class Impl;

private:
    const Impl& getImpl() const
    {
        return *m_Impl;
    }

    Impl& updImpl()
    {
        return *m_Impl;
    }

    std::shared_ptr<Impl> m_Impl = nullptr;

    friend CurveSegment;
    friend CurveSegment::Impl;
    friend CableSubsystem;

    // Befriend the helper class for testing the implementation.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                                CABLE SUBSYSTEM
//==============================================================================

class SimTK_SIMBODY_EXPORT CableSubsystem : public Subsystem
{
public:
    CableSubsystem();
    explicit CableSubsystem(MultibodySystem&);

    int getNumCables() const;
    const CableSpan& getCable(CableSpanIndex idx) const;
    CableSpan& updCable(CableSpanIndex idx);

    SimTK_PIMPL_DOWNCAST(CableSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;

    // Befriend the helper class for testing the implementation.
    friend CableSubsystemTestHelper;
};

//==============================================================================
//                      SUBSYSTEM TESTING HELPER
//==============================================================================

// Helper class for testing the internally computed path error jacboian.
class SimTK_SIMBODY_EXPORT CableSubsystemTestHelper
{
public:
    CableSubsystemTestHelper() = default;

    // Verify the computed path error jacobian by applying a small correction
    // to each CurveSegment's geodesic (i.e. a perturbation), and computing the
    // resulting change in the path error vector.
    bool applyPerturbationTest(
        const State& s,
        const CableSubsystem& subsystem,
        Real perturbation,
        Real bound,
        std::ostream& os);
};

} // namespace SimTK

#endif
