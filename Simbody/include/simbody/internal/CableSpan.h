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

SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanObstacleIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanIndex);

class MultibodySystem;
class CableSubsystem;
class CableSubsystemTestHelper;
class CableSpan;

//==============================================================================
//                          CABLE SPAN
//==============================================================================
// Representation of a cable spanning between two poins.
// The cable can adopt surface obstacles that must be wrapped over.
//
// NOTE: The interaction with the obstacles is **ordered**. The cable can wrap
// over the obstacles in the order that they are added. This means that the
// cable can not wrap over the first obstacle twice, even though it might
// spatially intersect it twice.
class SimTK_SIMBODY_EXPORT CableSpan final
{
public:
    using FrenetFrame   = Transform;
    using ObstacleIndex = CableSpanObstacleIndex;

    /** Identifies whether the cable is in contact with the obstacle surface or
     * not.
     * TODO the InitialGuess case is added as a fix for initializing the path,
     * and should be removed later. */
    enum class ObstacleWrappingStatus
    {
        InitialGuess, // Hotfix for path initialization.
        InContactWithSurface,
        LiftedFromSurface,
        Disabled, // Allow user to manually disable the obstacle.
    };

//------------------------------------------------------------------------------

    CableSpan();
    ~CableSpan();

    /** Clones the data, but not the subsystem entry. */
    CableSpan(const CableSpan& source);
    /** Clones the data, but not the subsystem entry. */
    CableSpan& operator=(const CableSpan& source);

    CableSpan(CableSpan&&) noexcept            = default;
    CableSpan& operator=(CableSpan&&) noexcept = default;

//------------------------------------------------------------------------------

    // Construct a new cable, and add it to the subsystem.
    // The cable spans from the origin point on the origin body, to the
    // terminal point on the terminal body.
    CableSpan(
        CableSubsystem& subsystem,
        MobilizedBodyIndex originBody,
        const Vec3& defaultOriginPoint, // TODO rename
        MobilizedBodyIndex terminationBody,
        const Vec3& defaultTerminationPoint);

//------------------------------------------------------------------------------
//                     OBSTACLE CONFIGURATION
//------------------------------------------------------------------------------

    /** Add an obstacle to the cable's path.
     * @param mobod The body that the contact geometry is rigidly attached to.
     * @param X_BS Transform specifying the location and orientation of the
     * contact geometry's origin frame with respect to the mobilized body.
     * @param geometry The contact geometry over which this segment wraps.
     * @param contactPointHint A guess of the contact point of the cable
     * span and the contact geometry to compute the initial cable path. This
     * point is defined in the local contact geometry's frame. The point will be
     * used as a starting point when computing the initial cable path. As such,
     * it does not have to lie on the contact geometry's surface, nor does it
     * have to belong to a valid cable path. */
    ObstacleIndex addSurfaceObstacle(
        MobilizedBodyIndex mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint);

    /** Get the number of obstacles added to the path. */
    int getNumSurfaceObstacles() const;

    // Helper function: TODO remove?
    const MobilizedBody& getObstacleMobilizedBody(ObstacleIndex ix) const;

    /** Get the index of the mobilized body that the obstacle is attached to. */
    const MobilizedBodyIndex& getObstacleMobilizedBodyIndex(
        ObstacleIndex ix) const;
    /** Set the index of the mobilized body that the obstacle is attached to. */
    void setObstacleMobilizedBodyIndex(
        ObstacleIndex ix,
        MobilizedBodyIndex body);

    /** Get the orientation and position of the obstacle's surface with respect
     * to its mobilized body. */
    const Transform& getObstacleXformSurfaceToBody(ObstacleIndex ix) const;
    /** Set the orientation and position of the obstacle's surface with respect
     * to its mobilized body. */
    void setObstacleXformSurfaceToBody(ObstacleIndex ix, Transform X_BS);

    /** Get the ContactGeometry attached to the obstacle */
    const ContactGeometry& getObstacleContactGeometry(ObstacleIndex ix) const;
    /** Set the ContactGeometry attached to the obstacle. TODO: dont take
     * ownership. */
    void setObstacleContactGeometry(ObstacleIndex ix, ContactGeometry geometry);

    /** Get the point on the obstacle used to compute the initial path. */
    Vec3 getObstacleInitialContactPointHint(ObstacleIndex ix) const;
    /** Set the point on the obstacle used to compute the initial path. */
    void setObstacleInitialContactPointHint(
        ObstacleIndex ix,
        Vec3 initialContactPointHint);

//------------------------------------------------------------------------------
//                     OBSTACLE CALCULATIONS
//------------------------------------------------------------------------------

    /** Get the wrapping status of the cable path over the given obstacle.
     * State must be realized to Stage::Position. */
    ObstacleWrappingStatus getObstacleWrappingStatus(
        const State& state,
        ObstacleIndex ix) const;

    /** Get the length of the curve segment over the given obstacle.
     * Throws an exception if the cable is not in contact with the surface.
     * State must be realized to Stage::Position. */
    Real getCurveSegmentLength(const State& state, ObstacleIndex ix) const;

    int getCurveSegmentNumberOfIntegratorStepsTaken(
        const State& state,
        ObstacleIndex ix) const;
    Real getCurveSegmentInitialIntegratorStepSize(
        const State& state,
        ObstacleIndex ix) const;

    /** Compute points along the obstacle's curve segment at equal length
     * intervals.
     * The computed points are interpolated from the underlying geodesic, and
     * might slightly violate configured tolerances.
     * State must be realized to Stage::Position.
     *
     * @param state System State.
     * @param ix Index of the obstacle that the curve wraps over.
     * @param nPoints Number of points used to resample the curve. Minimum
     * allowed is two.
     * @param sink Where the interpolated point will be written to. Both the
     * length measured from the curve segment's first contact point, as well as
     * the position are written.
     * @return Number of points actually written. If the curve has zero length,
     * only one point will be computed. If the cable is not in contact with
     * this obstacle, zero points will be written. */
    int calcCurveSegmentPathPoints(
        const State& state,
        ObstacleIndex ix,
        int nPoints,
        std::function<void(Real length, Vec3 point)> sink) const;

    // TODO same as calcCurveSegmentPathPoints but also computes the curve's
    // tangent.
    int calcCurveSegmentPathPointsAndTangents(
        const State& state,
        ObstacleIndex ix,
        int nPoints,
        std::function<void(Real length, Vec3 point, UnitVec3 tangent)> sink)
        const;

//------------------------------------------------------------------------------
//                     CABLE CONFIGURATION
//------------------------------------------------------------------------------

    /** Get the tolerance used to enfore the surface constraints of all
     * obstacles. */
    Real getSurfaceConstraintTolerance() const;
    /** Set the tolerance used to enfore the surface constraints of all
     * obstacles. */
    void setSurfaceConstraintTolerance(Real tolerance);

    /** Get the maximum number of solver iterations allowed when enforcing the
     * surface constraints of all obstacles. */
    int getSurfaceProjectionMaxIter()
        const; // TODO not connected to anything currently.
    /** Set the maximum number of solver iterations allowed when enforcing the
     * surface constraints of all obstacles. */
    void setSurfaceProjectionMaxIter(int maxIter);

    /** Get the accuracy used by the numerical integrator when computing a
     * geodesic over an obstacle. */
    Real getIntegratorAccuracy() const;
    /** Set the accuracy used by the numerical integrator when computing a
     * geodesic over an obstacle. */
    void setIntegratorAccuracy(Real accuracy);

    /** Get the maximum number of solver iterations for finding the optimal
     * path. */
    int getSolverMaxIter() const;
    /** Set the maximum number of solver iterations for finding the optimal
     * path. */
    void setSolverMaxIter(int maxIter);

    // TODO merge with surface constraint tolerance above?
    Real getPathErrorAccuracy() const;
    void setPathErrorAccuracy(Real accuracy);

    // TODO advanced parameter: Max allowed step during Newton iteration in
    // angles. Using the local curvatures the max angular step is converted to a
    // max allowed translation.
    Real getMaxRadialStepInDegrees() const;
    void setMaxRadialStepInDegrees(Real maxStepInDeg);

//------------------------------------------------------------------------------
//                     CABLE CALCULATIONS
//------------------------------------------------------------------------------

    /** Get the total cable length.
     * State must be realized to Stage::Position. */
    Real getLength(const State& state) const;

    /** Get the derivative of the total cable length.
     * State must be realized to Stage::Position. */
    Real getLengthDot(const State& state) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    /** Compute the cable power.
     * State must be realized to Stage::Position. */
    Real calcCablePower(const State& state, Real tension) const;

    /** Compute points on the path spanned by this cable.
     * The curve segments will be resampled from the underlying geodesic such
     * that the computed points are at equal length intervals. The straight
     * line segments are not resampled.
     * State must be realized to Stage::Position.
     *
     * @param state System State.
     * @param lengthInterval Each curve segments will be resampled such the the
     * points interspacing does not exceed this value.
     * @param sink Where the path points will be written to.
     * @return The total number of computed points.
     */
    int calcPathPoints(
        const State& state,
        Real lengthInterval,
        std::function<void(Real length, Vec3 point)> sink) const;

    /** Number of solver iterations required to compute the cable's path.
    * State must be realized to Stage::Position. */
    int getNumSolverIter(const State& state) const;

    /** Maximum path error of the current cable's path.
    * TODO explain path error.
    * State must be realized to Stage::Position. */
    Real getMaxPathError(const State& state) const;

    /** TODO Remove this?
    * Overwrite the path used as a warmstart for the solver. This path is
    * normally auto-updated after completing an integrator step. */
    void storeCurrentPath(State& state) const;

    class Impl;

private:
    /** Wrap the Impl, increasing the reference count. */
    explicit CableSpan(std::shared_ptr<Impl> impl): m_Impl(std::move(impl)) {}

    /** Cheap copy of pointer to the Impl. */
    CableSpan copyImpl() const {
        return CableSpan(m_Impl);
    }

    const Impl& getImpl() const
    {
        return *m_Impl;
    }

    Impl& updImpl()
    {
        return *m_Impl;
    }

    std::shared_ptr<Impl> m_Impl = nullptr;

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
