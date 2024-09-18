#ifndef SimTK_SIMBODY_CABLE_SPAN_H_
#define SimTK_SIMBODY_CABLE_SPAN_H_

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

#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/common.h"

namespace SimTK
{

/** @class SimTK::CableSpanObstacleIndex
This is a unique integer type for identifying obstacles comprising a particular
CableSpan. These begin at zero for each CableSpan. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanObstacleIndex);

/** @class SimTK::CableSpanIndex
This is a unique integer type for quickly identifying specific cables for fast
lookup purposes. These begin at zero for each CableSubsystem. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanIndex);

class MultibodySystem;
class CableSubsystem;
class CableSubsystemTestHelper;

//==============================================================================
/** This class represents the path of a frictionless cable from an origin point
fixed to a body, over geometric obstacles fixed to other bodies, to a final
termination point.

The CableSpan's path can be seen as consisting of straight line segments and
curved line segments: A curved segment over each obstacle, and straight segments
connecting them to each other and to the end points. Each curved segment is
computed as a geodesic to give (in some sense) a shortest path over the surface.
During a simulation the cable can slide freely over the obstacle surfaces. It
can lose contact with a surface and miss that obstacle in a straight line.
Similarly the cable can touchdown on the obstacle if the surface obstructs the
straight line segment again.

The path is computed as an optimization problem using the previous optimal path
as the warm start. This is done by computing natural geodesic corrections for
each curve segment to compute the locally shortest path, as described in the
following publication:

    Scholz, A., Sherman, M., Stavness, I. et al (2015). A fast multi-obstacle
    muscle wrapping method using natural geodesic variations. Multibody System
    Dynamics 36, 195–219.

The overall path is locally the shortest, allowing winding over an obstacle
multiple times, without flipping to the other side.

During initialization the path is assumed to be in contact with each obstacle at
the user defined contact-point-hint. At each contact point a
zero-length-geodesic is computed with the tangent estimated as pointing from the
previous path point to the next path point. From this configuration, the solver
is started, and the geodesics will be corrected until the entire path is smooth.

Important to note is that the cable's interaction with the obstacles is ordered
based on the order in which they were added. That is, if three obstacles are
added to a cable, then, the cable can wrap over the first, then the second, and
then the third. If the first obstacle spatially collides with the path twice it
will not actually wrap over it twice. Use CablePath if this is not the desired
behavior.

Note that a CableSpan is a geometric object, not a force or constraint element.
That is, a CableSpan alone will not influence the behavior of a simulation.
However, forces and constraint elements can be constructed that make use of a
CableSpan to generate forces.

A CableSpan must be registered with a CableSubsystem which manages their
runtime evaluation. **/
class SimTK_SIMBODY_EXPORT CableSpan final {
public:
    /** Construct a new cable that can be configured later. **/
    CableSpan();
    /** Delete the cable if this handle was the last reference to it. **/
    ~CableSpan() noexcept = default;
    /** Copy constructor is shallow and reference counted. **/
    CableSpan(const CableSpan&) = default;
    /** Copy assignment is shallow and reference counted. **/
    CableSpan& operator=(const CableSpan& source) = default;
    CableSpan(CableSpan&&) noexcept               = default;
    CableSpan& operator=(CableSpan&&) noexcept    = default;

    /** @name Cable construction */
    ///@{

    /** Create a straight-line cable spanning between a point fixed on one body
    to one fixed on another body. You can add obstacles and move the end points
    later.
    @param subsystem The subsystem that this cable is adopted by.
    @param originBody The mobilized body that the origin point is rigidly
    attached to.
    @param originStation The origin point of the cable defined in body fixed
    coordinates.
    @param terminationBody The mobilized body that the termination point is
    rigidly attached to.
    @param terminationStation The termination point of the cable defined in
    body fixed coordinates. **/
    CableSpan(
        CableSubsystem& subsystem,
        MobilizedBodyIndex originBody,
        const Vec3& originStation,
        MobilizedBodyIndex terminationBody,
        const Vec3& terminationStation);

    /** Add an obstacle to the cable's path that must be wrapped over.
    @param obstacleBody The body that the contact geometry is rigidly attached
    to.
    @param X_BS Transform specifying the location and orientation of the
    contact geometry's origin frame with respect to the mobilized body.
    @param obstacleGeometry The geometry of the obstacle's surface.
    @return The index of the added obstacle in this cable. **/
    CableSpanObstacleIndex addObstacle(
        MobilizedBodyIndex obstacleBody,
        const Transform& X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry);

    /** Add an obstacle to the cable's path that must be wrapped over.
    @param obstacleBody The body that the contact geometry is rigidly attached
    to.
    @param X_BS Transform specifying the location and orientation of the
    contact geometry's origin frame with respect to the mobilized body.
    @param obstacleGeometry The geometry of the obstacle's surface.
    @param contactPointHint_S A guess for the cable's initial contact point on
    the obstacle, in local surface frame coordinates. The point will be used as
    a starting point when computing the initial cable path. As such, it does
    not have to lie on the contact geometry's surface, nor does it have to
    belong to a valid cable path.
    @return The index of the added obstacle in this cable. **/
    CableSpanObstacleIndex addObstacle(
        MobilizedBodyIndex obstacleBody,
        const Transform& X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry,
        const Vec3& contactPointHint_S);

    /** Get the number of obstacles added to the path. **/
    int getNumObstacles() const;

    ///@}

    /** @name End points configuration */
    ///@{

    /** Get the index of the mobilized body that the cable's origin point is
    attached to. **/
    MobilizedBodyIndex getOriginBodyIndex() const;

    /** Set the index of the mobilized body that the cable's origin point is
    attached to. **/
    void setOriginBodyIndex(MobilizedBodyIndex originBody);

    /** Get the index of the mobilized body that the cable's termination point
    is attached to. **/
    MobilizedBodyIndex getTerminationBodyIndex() const;

    /** Set the index of the mobilized body that the cable's termination point
    is attached to. **/
    void setTerminationBodyIndex(MobilizedBodyIndex terminationBody);

    /** Get the cable's origin point defined in body fixed coordinates. **/
    Vec3 getOriginStation() const;

    /** Set the cable's origin point defined in body fixed coordinates. **/
    void setOriginStation(const Vec3& originStation);

    /** Get the cable's termination point defined in body fixed coordinates.
    **/
    Vec3 getTerminationStation() const;

    /** Set the cable's termination point defined in body fixed coordinates.
    **/
    void setTerminationStation(const Vec3& terminationStation);

    ///@}

    /** @name Obstacle configuration */
    ///@{

    /** Get the index of the mobilized body that the obstacle is attached to.
    @param ix The index of the obstacle in this CableSpan.
    @return The index of the mobilized body that the obstacle is attached
    to. **/
    const MobilizedBodyIndex& getObstacleMobilizedBodyIndex(
        CableSpanObstacleIndex ix) const;

    /** Set the index of the mobilized body that the obstacle is attached to.
    @param ix The index of the obstacle in this CableSpan.
    @param body The index of the mobilized body that the obstacle is attached
    to. **/
    void setObstacleMobilizedBodyIndex(
        CableSpanObstacleIndex ix,
        MobilizedBodyIndex body);

    /** Get the orientation and position of the obstacle's surface with respect
    to its mobilized body.
    @param ix The index of the obstacle in this CableSpan.
    @return The orientation and position of the obstacle's surface with respect
    to it's mobilized body. **/
    const Transform& getObstacleTransformSurfaceToBody(
        CableSpanObstacleIndex ix) const;

    /** Set the orientation and position of the obstacle's surface with respect
    to its mobilized body.
    @param ix The index of the obstacle in this CableSpan.
    @param X_BS The orientation and position of the obstacle's surface with
    respect to it's mobilized body. **/
    void setObstacleTransformSurfaceToBody(
        CableSpanObstacleIndex ix,
        const Transform& X_BS);

    /** Get the obstacle's ContactGeometry.
    @param ix The index of the obstacle in this CableSpan.
    @return The obstacle's surface geometry. **/
    const ContactGeometry& getObstacleContactGeometry(
        CableSpanObstacleIndex ix) const;

    /** Set the obstacle's ContactGeometry.
    @param ix The index of the obstacle in this CableSpan.
    @param obstacleGeometry The obstacle's surface geometry. **/
    void setObstacleContactGeometry(
        CableSpanObstacleIndex ix,
        std::shared_ptr<const ContactGeometry> obstacleGeometry);

    /** Get the point on the obstacle used to compute the initial path.
    @param ix The index of the obstacle in this CableSpan.
    @return The point in the obstacle's local surface coordinates. **/
    Vec3 getObstacleContactPointHint(CableSpanObstacleIndex ix) const;

    /** Set the point on the obstacle used to compute the initial path.
    @param ix The index of the obstacle in this CableSpan.
    @param contactPointHint_S A guess for the cable's initial contact point on
    the obstacle, in local surface frame coordinates. The point will be used as
    a starting point when computing the initial cable path. As such, it does
    not have to lie on the contact geometry's surface, nor does it have to
    belong to the optimal cable path. **/
    void setObstacleContactPointHint(
        CableSpanObstacleIndex ix,
        Vec3 contactPointHint_S);

    ///@}

    /** @name Solver configuration */
    ///@{

    /** Get the accuracy used by the numerical integrator when computing a
    geodesic over an obstacle.
    Note: This does not affect the integrator that is used to propagate the
    multibody system over time, that is a different integrator. **/
    Real getCurveSegmentAccuracy() const;

    /** Set the accuracy used by the numerical integrator when computing a
    geodesic over an obstacle.
    Note: This does not affect the integrator that is used to propagate the
    multibody system over time, that is a different integrator. **/
    void setCurveSegmentAccuracy(Real accuracy);

    /** Get the maximum number of solver iterations for finding the optimal
    path. **/
    int getSolverMaxIterations() const;

    /** Set the maximum number of solver iterations for finding the optimal
    path. **/
    void setSolverMaxIterations(int maxIterations);

    /** Number of solver iterations used to compute the current cable's path.
    State must be realized to Stage::Position.
    @param state System State.
    @return Number of solver iterations used. **/
    int getNumSolverIterations(const State& state) const;

    /** Get the smoothness tolerance used to compute the optimal path.
    The (non) smoothness is defined as the angular discontinuity at the points
    where straight- and curved-line segments meet, measured in radians. When
    computing the optimal path this smoothness is optimized, and the solver
    stops when reaching given tolerance.
    **/
    Real getSmoothnessTolerance() const;

    /** Set the smoothness tolerance used to compute the optimal path.
    See CableSpan::getSmoothnessTolerance.
    **/
    void setSmoothnessTolerance(Real tolerance);

    /** The smoothness of the current cable's path.
    See CableSpan::getSmoothnessTolerance.
    State must be realized to Stage::Position. **/
    Real getSmoothness(const State& state) const;

    ///@}

    /** @name Path computations */
    ///@{

    /** Get the total cable length.
    State must be realized to Stage::Position.
    @param state State of the system.
    @return The total cable length. **/
    Real calcLength(const State& state) const;

    /** Get the derivative of the total cable length.
    State must be realized to Stage::Velocity.
    @param state State of the system.
    @return The time derivative of the total cable length. **/
    Real calcLengthDot(const State& state) const;

    /** Given a tension > 0 acting uniformly along this cable, compute the
    resulting forces applied to the bodies it touches.
    The body forces are added into the appropriate slots in the supplied Array
    which has one entry per body in the same format as is supplied to the
    calcForce() method of force elements. If the supplied tension is <= 0,
    signifying a slack cable, this method does nothing. **/
    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    /** Compute the cable power.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param tension The tension force in the cable.
    @return The cable power. **/
    Real calcCablePower(const State& state, Real tension) const;

    /** Compute points on the path spanned by this cable suitable for
    visualization purposes.
    State must be realized to Stage::Position.
    @param state System State.
    @param sink Where the path points (in ground frame) will be written to. **/
    void calcDecorativePathPoints(
        const State& state,
        const std::function<void(Vec3 point_G)>& sink) const;

    ///@}

    /** @name Curve segment computations */
    ///@{

    /** Returns true when the cable is in contact with the obstacle.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan. **/
    bool isInContactWithObstacle(const State& state, CableSpanObstacleIndex ix)
        const;

    /** Compute the Frenet frame associated with the obstacle's curve segment at
    the initial contact point on that obstacle.
    If the path is not in contact with the obstacle's surface the frame will
    contain invalid data (NaNs). The Frenet frame is measured relative to
    ground, with the tangent along the X axis, the surface normal along the Y
    axis and the binormal along the Z axis.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan.
    @return The Frenet frame at the obstacle's initial contact point. **/
    Transform calcCurveSegmentInitialFrenetFrame(
        const State& state,
        CableSpanObstacleIndex ix) const;

    /** Compute the Frenet frame associated with the obstacle's curve segment at
    the final contact point on that obstacle.
    If the path is not in contact with the obstacle's surface the frame will
    contain invalid data (NaNs). The Frenet frame is measured relative to
    ground, with the tangent along the X axis, the surface normal along the Y
    axis and the binormal along the Z axis.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan.
    @return The Frenet frame at the obstacle's final contact point. **/
    Transform calcCurveSegmentFinalFrenetFrame(
        const State& state,
        CableSpanObstacleIndex ix) const;

    /** Get the arc length of the obstacle's curve segment.
    Returns NaN if the obstacle is not in contact with the path.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan.
    @return The arc length. **/
    Real calcCurveSegmentArcLength(
        const State& state,
        CableSpanObstacleIndex ix) const;

    /** Compute points along the obstacle's curve segment at equal length
    intervals.

    The first and last contact point of the curve segment are always included.
    Therefore the minimum allowed number of samples is 2.

    For analytic surfaces the resampled points are computed analytically.
    Otherwise, the points are computed from the internally stored geodesic
    using Hermite interpolation, generally resulting in a good approximation of
    the point at the requested arc lengths.

    This function throws an exception if the obstacle is not in contact with
    the cable. Check isInContactWithObstacle() before calling this function.

    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan.
    @param numSamples The number of resampling points (minimum is 2).
    @param sink Where the resampled points (in ground frame) will be written
    to. **/
    void calcCurveSegmentResampledPoints(
        const State& state,
        CableSpanObstacleIndex ix,
        int numSamples,
        const std::function<void(Vec3 point_G)>& sink) const;

    ///@}

    /** @cond **/ // Hide from Doxygen.
    class Impl;

private:
    const Impl& getImpl() const
    {
        return *m_impl;
    }

    Impl& updImpl()
    {
        return *m_impl;
    }

    std::shared_ptr<Impl> m_impl;

    friend CableSubsystem;

    // Befriend the helper class for testing the implementation.
    friend CableSubsystemTestHelper;
    /** @endcond **/
};

//==============================================================================
/** This subsystem manages cables spanning between two points in a system.
Each cable is represented by a CableSpan. Obstacles can be added to the cable,
and must be wrapped over. The calculated path will consist of a series of
straight line segments between obstacles, and curved segments (geodesics) over
the obstacles. Force elements defined elsewhere may make use of the cables to
apply forces to the system.

See CableSpan. **/
class SimTK_SIMBODY_EXPORT CableSubsystem : public Subsystem {
public:
    explicit CableSubsystem(MultibodySystem&);

    const MultibodySystem& getMultibodySystem() const;

    /** Get the number of cables being managed by this subsystem.
    These are identified by CableSpanIndex values from 0 to getNumCables()-1.
    This is available after realizeTopology() and does not change subsequently.
    **/
    int getNumCables() const;

    /** Get const access to a particular cable. **/
    const CableSpan& getCable(CableSpanIndex ix) const;

    /** Get writeable access to a particular cable. **/
    CableSpan& updCable(CableSpanIndex ix);

    /** @cond **/ // Hide from Doxygen.
    SimTK_PIMPL_DOWNCAST(CableSubsystem, Subsystem);

    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;

    // Befriend the helper class for testing the implementation.
    friend CableSubsystemTestHelper;
    /** @endcond **/
};

/** @cond **/ // Hide from Doxygen.
//==============================================================================
/** Helper class for testing correctness of all currently computed paths in a
CableSubsystem.

For each CableSpan's computed path it is checked that:
1. Any curve segment over an obstacle is a correct geodesic.
2. The path error Jacobian matches a perturbation test.

The geodesics and Jacobians are internal to the CableSpan's implementation, we
refer to the following publication to understand their role:

    Scholz, A., Sherman, M., Stavness, I. et al (2015). A fast multi-obstacle
    muscle wrapping method using natural geodesic variations. Multibody System
    Dynamics 36, 195–219.

Consider a bug in the code, then this test is designed to have
high sensitivity (and lower specificity). That is, if you pass, you can rest
assured that all curve segments are indeed geodesics, and that the Jacobian
correctly predicts the local effect of the NaturalGeodesicCorrection on the
path error vector. It is possible to fail this test in the absence of any bugs
by using incompatible configuration parameters; e.g. by setting a very poor
integrator accuracy for the CableSpan, and using a very small perturbation value
in the Jacobian perturbation test. For use in a unit test a high sensitivity
works just fine: If you pass it normally, you should pass it after a code
refactor as well.

Testing each computed path in your simulation is possible, but will criple your
simulation speed.

State must be realized to Stage::Position before testing. **/
class SimTK_SIMBODY_EXPORT CableSubsystemTestHelper {
public:
    /** Construct a CableSubsystemTestHelper that can be used to test a
    CableSubsystem. **/
    CableSubsystemTestHelper();
    ~CableSubsystemTestHelper() noexcept;
    CableSubsystemTestHelper(const CableSubsystemTestHelper&);
    CableSubsystemTestHelper& operator=(const CableSubsystemTestHelper&);
    CableSubsystemTestHelper(CableSubsystemTestHelper&&) noexcept = default;
    CableSubsystemTestHelper& operator=(CableSubsystemTestHelper&&) noexcept =
        default;

    /** Test the paths of all CableSpan registered in a CableSubsystem.
    State must be realized to Stage::Position.
    An exception is thrown when the test fails.
    @param state State of the system.
    @param subsystem CableSubsystem containing all paths to be tested.
    @param testReport Interesting test results will be written to this stream.
    **/
    void testCurrentPath(
        const State& state,
        const CableSubsystem& subsystem,
        std::ostream& testReport) const;

    class Impl;

private:
    std::unique_ptr<Impl> m_impl;
};
/** @endcond **/

} // namespace SimTK

#endif
