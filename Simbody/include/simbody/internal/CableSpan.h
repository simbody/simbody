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
    using FrenetFrame = Transform;

    enum class ObstacleWrappingStatus {
        InitialGuess,
        InContactWithSurface,
        LiftedFromSurface,
        Disabled,
    };

//------------------------------------------------------------------------------

    // TODO who owns the cable span? the subsystem?
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
//                     OBSTACLE CONFIGURATION
//------------------------------------------------------------------------------

    CurveSegmentIndex addSurfaceObstacle(
        MobilizedBodyIndex mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint);

    int getNumSurfaceObstacles() const;

    const MobilizedBody& getObstacleMobilizedBody(CurveSegmentIndex ix) const;
    const MobilizedBodyIndex& getObstacleMobilizedBodyIndex(CurveSegmentIndex ix) const;
    void setObstacleMobilizedBodyIndex(CurveSegmentIndex ix, MobilizedBodyIndex body);

    const Transform& getObstacleXformSurfaceToBody(CurveSegmentIndex ix) const;
    void setObstacleXformSurfaceToBody(CurveSegmentIndex ix, Transform X_BS);

    const ContactGeometry& getObstacleContactGeometry(CurveSegmentIndex ix) const;
    void setObstacleContactGeometry(CurveSegmentIndex ix, ContactGeometry geometry);

    Vec3 getObstacleInitialContactPointHint(CurveSegmentIndex ix) const;
    void setObstacleInitialContactPointHint(CurveSegmentIndex ix, Vec3 initialContactPointHint);

//------------------------------------------------------------------------------
//                     OBSTACLE CALCULATIONS
//------------------------------------------------------------------------------

    ObstacleWrappingStatus getObstacleWrappingStatus(const State& state, CurveSegmentIndex ix) const;
    Real getCurveSegmentLength(const State& state, CurveSegmentIndex ix) const;

    int calcCurveSegmentPathPoints(const State& state, CurveSegmentIndex ix, int nPoints, std::function<void(Vec3 point)>& sink) const;

    // This is useful for debugging and visualization.
    int calcCurveSegmentPathPointsAndTangents(const State& state, CurveSegmentIndex ix, int nPoints, std::function<void(Vec3 point, UnitVec3 tangent)>& sink) const;

    int getCurveSegmentNumberOfIntegratorStepsTaken(const State& state, CurveSegmentIndex ix) const;
    Real getCurveSegmentInitialIntegratorStepSize(const State& state, CurveSegmentIndex ix) const;

    // TODO perhaps remove this, and replace with calcTangentAtLength(), calcPosAtLength().
    const FrenetFrame& getCurveSegmentFirstFrenetFrame(const State& state, CurveSegmentIndex ix) const;
    // TODO perhaps remove this.
    const FrenetFrame& getCurveSegmentLastFrenetFrame(const State& state, CurveSegmentIndex ix) const;

//------------------------------------------------------------------------------
//                     CABLE CONFIGURATION
//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------
//                     CABLE CALCULATIONS
//------------------------------------------------------------------------------

    Real getLength(const State& state) const;

    Real getLengthDot(const State& state) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    Real calcCablePower(const State& state, Real tension) const;

    int calcPathPoints(const State& state, Real maxLengthIncrement, std::function<void(Vec3 point)>& sink) const;

    class Impl;

private:
    explicit CableSpan(Impl);

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
