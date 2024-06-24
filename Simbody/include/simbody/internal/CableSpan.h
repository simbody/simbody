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
CableSpan. These begin at zero for each cable. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanObstacleIndex);

/** @class SimTK::CableSpanIndex
This is a unique integer type for quickly identifying specific cables for fast
lookup purposes. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanIndex);

class MultibodySystem;
class CableSubsystem;
class CableSubsystemTestHelper;

//==============================================================================
//                                CABLESPAN
//==============================================================================
/** Class CableSpan models a cable spanning between two points with obstacles
that must be wrapped over.

TODO add description (similar to CablePath, depending on how things develop).

NOTE: The interaction with the obstacles is **ordered**. The cable can wrap over
the obstacles in the order that they are added. This means that the cable can
not wrap over the first obstacle twice, even though it might spatially intersect
it twice.

Algorithm details can be found in the following publication:

    Scholz, A., Sherman, M., Stavness, I. et al (2016). A fast multi-obstacle
    muscle wrapping method using natural geodesic variations. Multibody System
    Dynamics 36, 195–219.

**/
class SimTK_SIMBODY_EXPORT CableSpan final {
public:
    /** Construct a new cable that can be configured later. **/
    CableSpan();
    /** Delete the cable if this handle was the last reference to it. **/
    ~CableSpan() noexcept = default;
    /** Copy constructor is shallow and reference counted. **/
    CableSpan(const CableSpan& source) = default;
    /** Copy assignment is shallow and reference counted. **/
    CableSpan& operator=(const CableSpan& source) = default;
    CableSpan(CableSpan&&) noexcept               = default;
    CableSpan& operator=(CableSpan&&) noexcept    = default;

    /** Create a straight-line cable spanning between a point fixed on one body
    to one fixed on another body. You can add additional obstacles and move the
    end points later.
    @param subsystem The subsystem that this cable is adopted by.
    @param originBody The mobilized body that the origin point is rigidly
    attached to.
    @param originPoint_B The origin point of the cable defined in body fixed
    coordinates.
    @param terminationBody The mobilized body that the termination point is
    rigidly attached to.
    @param terminationPoint_B The termination point of the cable defined in
    body fixed coordinates. **/
    CableSpan(
        CableSubsystem& subsystem,
        MobilizedBodyIndex originBody,
        const Vec3& originPoint_B,
        MobilizedBodyIndex terminationBody,
        const Vec3& terminationPoint_B);

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
    CableSpanObstacleIndex addSurfaceObstacle(
        MobilizedBodyIndex obstacleBody,
        const Transform& X_BS,
        std::shared_ptr<const ContactGeometry> obstacleGeometry,
        Vec3 contactPointHint_S);

    /** Get the number of obstacles added to the path. **/
    int getNumSurfaceObstacles() const;

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
    const Transform& getObstacleXformSurfaceToBody(
        CableSpanObstacleIndex ix) const;

    /** Set the orientation and position of the obstacle's surface with respect
    to its mobilized body.
    @param ix The index of the obstacle in this CableSpan.
    @param X_BS The orientation and position of the obstacle's surface with
    respect to it's mobilized body. **/
    void setObstacleXformSurfaceToBody(
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
    Vec3 getObstacleInitialContactPointHint(CableSpanObstacleIndex ix) const;

    /** Set the point on the obstacle used to compute the initial path.
    @param ix The index of the obstacle in this CableSpan.
    @param contactPointHint_S A guess for the cable's initial contact point on
    the obstacle, in local surface frame coordinates. The point will be used as
    a starting point when computing the initial cable path. As such, it does
    not have to lie on the contact geometry's surface, nor does it have to
    belong to a valid cable path. **/
    void setObstacleInitialContactPointHint(
        CableSpanObstacleIndex ix,
        Vec3 contactPointHint_S);

    /** Returns true when the cable is in contact with the obstacle.
    State must be realized to Stage::Position.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan. **/
    bool isInContactWithObstacle(const State& state, CableSpanObstacleIndex ix)
        const;

    /** Compute knot points of the curve segment (i.e. the geodesic) on the
    obstacle in local coordinates, together with the transform for conversion to
    ground frame. If the cable is not in contact with the given obstacle,
    nothing will be written.

    TODO I'm unsure what kind of info we would like from the geodesics.
    TODO it feels like access to some info on the geodesic would be nice.
    TODO Feedback is welcome!
    **/
    void calcCurveSegmentKnots(
        const State& state,
        CableSpanObstacleIndex ix,
        const std::function<void(
            const ContactGeometry::GeodesicKnotPoint& geodesicKnot_S,
            const Transform& X_GS)>& sink) const;

    /** Disable contact between the obstacle and this CableSpan.
    State must be realized to Stage::Instance.
    Does nothing if the obstacle was already disabled. Otherwise this method
    will invalidate Stage::Position and later stages.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan. **/
    void setObstacleContactDisabled(
        const State& state,
        CableSpanObstacleIndex ix) const;

    /** Enable contact between the obstacle and this CableSpan.
    State must be realized to Stage::Instance.
    Does nothing if the obstacle was already enabled. Otherwise this method
    will invalidate Stage::Position and later stages.
    @param state State of the system.
    @param ix The index of the obstacle in this CableSpan. **/
    void setObstacleContactEnabled(
        const State& state,
        CableSpanObstacleIndex ix) const;

    /** Get the tolerance used to enfore the surface constraints of all
    obstacles. **/
    Real getSurfaceConstraintTolerance() const;

    /** Set the tolerance used to enfore the surface constraints of all
    obstacles.

    TODO Currently the constraint projection tolerance is set for the entire
    path. Which seems very convenient. However, this assumes that the surface
    constraints of all obstacles are defined on the same scale. Which is
    generally not true. Perhaps it is worth defining the surface constraints on
    a similar scale. **/
    void setSurfaceConstraintTolerance(Real tolerance);

    /** Get the accuracy used by the numerical integrator when computing a
    geodesic over an obstacle. **/
    Real getIntegratorAccuracy() const;

    /** Set the accuracy used by the numerical integrator when computing a
    geodesic over an obstacle. **/
    void setIntegratorAccuracy(Real accuracy);

    /** Get the maximum number of solver iterations for finding the optimal
    path. **/
    int getSolverMaxIterations() const;

    /** Set the maximum number of solver iterations for finding the optimal
    path. **/
    void setSolverMaxIterations(int maxIterations);

    // TODO merge with surface constraint tolerance above?
    Real getPathErrorAccuracy() const;
    void setPathErrorAccuracy(Real accuracy);

    // TODO advanced parameter: Max allowed step during Newton iteration in
    // angles. Using the local curvatures the max angular step is converted to a
    // max allowed translation.
    Real getMaxGeodesicCorrectionInDegrees() const;
    void setMaxGeodesicCorrectionInDegrees(Real maxCorrectionInDegrees);

    /** Get the total cable length.
    State must be realized to Stage::Position.
    @param state State of the system.
    @return The total cable length. **/
    Real getLength(const State& state) const;

    /** Get the derivative of the total cable length.
    State must be realized to Stage::Velocity.
    @param state State of the system.
    @return The time derivative of the total cable length. **/
    Real getLengthDot(const State& state) const;

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

    /** Number of solver iterations used to compute the current cable's path.
    State must be realized to Stage::Position.
    @param state System State.
    @return Number of solver iterations used. **/
    int getNumSolverIter(const State& state) const;

    /** Maximum path error of the current cable's path.
     * TODO explain path error.
     * State must be realized to Stage::Position. **/
    Real getMaxPathError(const State& state) const;

    /** TODO Remove this?
     * Overwrite the path used as a warmstart for the solver. This path is
     * normally auto-updated after completing an integrator step. **/
    void storeCurrentPath(State& state) const;

    /** @cond **/ // Hide from Doxygen.
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

    std::shared_ptr<Impl> m_Impl;

    friend CableSubsystem;

    // Befriend the helper class for testing the implementation.
    friend CableSubsystemTestHelper;
    /** @endcond **/
};

//==============================================================================
//                                CABLE SUBSYSTEM
//==============================================================================
/** TODO add description (similar to CableTrackerSubsystem, depending on how
things develop). **/
class SimTK_SIMBODY_EXPORT CableSubsystem : public Subsystem {
public:
    CableSubsystem();
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

//==============================================================================
//                      SUBSYSTEM TESTING HELPER
//==============================================================================
// TODO move this elsewhere?
// Helper class for testing the internally computed path error jacboian.
class SimTK_SIMBODY_EXPORT CableSubsystemTestHelper {
public:
    CableSubsystemTestHelper() = default;

    // Verify the computed path error jacobian by applying a small correction
    // to each CurveSegment's geodesic (i.e. a perturbation), and computing the
    // resulting change in the path error vector.
    bool applyPerturbationTest(
        const MultibodySystem& system,
        const CableSubsystem& subsystem,
        const State& s,
        Real perturbation,
        Real bound,
        std::ostream& os);
};

} // namespace SimTK

#endif
