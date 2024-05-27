#ifndef SimTK_SIMBODY_CABLE_SPAN_IMPL_H_
#define SimTK_SIMBODY_CABLE_SPAN_IMPL_H_

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

#include "CableSpan_CurveSegment_Impl.h"
#include "simbody/internal/CableSpan.h"

namespace SimTK
{

//==============================================================================
//                         CABLESPAN :: IMPL
//==============================================================================

class CableSpan::Impl
{
public:
    // Position level cache entry.
    struct PosInfo final
    {
        // Positions of the cable end points w.r.t. Ground.
        Vec3 originPoint_G{NaN, NaN, NaN};
        Vec3 terminationPoint_G{NaN, NaN, NaN};

        // Tangents at the cable end points w.r.t. Ground.
        UnitVec3 originTangent_G{NaN, NaN, NaN};
        UnitVec3 terminationTangent_G{NaN, NaN, NaN};

        Real cableLength = NaN;

        // Number of iterations required by solver to compute the path.
        size_t loopIter = 0;
    };

    // Velocity level cache entry.
    struct VelInfo final
    {
        Real lengthDot = NaN;
    };

//------------------------------------------------------------------------------

    Impl(
        CableSubsystem& subsystem,
        MobilizedBodyIndex originBody,
        Vec3 originPoint,
        MobilizedBodyIndex terminationBody,
        Vec3 terminationPoint) :
        m_Subsystem(&subsystem),
        m_OriginBody(originBody), m_OriginPoint(originPoint),
        m_TerminationBody(terminationBody), m_TerminationPoint(terminationPoint)
    {}

    CurveSegmentIndex adoptSegment(const CurveSegment& segment)
    {
        m_CurveSegments.push_back(segment);
        return CurveSegmentIndex(m_CurveSegments.size() - 1);
    }

//------------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    void realizePosition(const State& state) const;
    void realizeVelocity(const State& state) const;
    void invalidateTopology()
    {
        getSubsystem().invalidateSubsystemTopologyCache();
    }
    void invalidatePositionLevelCache(const State& state) const;

    const PosInfo& getPosInfo(const State& state) const;
    const VelInfo& getVelInfo(const State& state) const;

//------------------------------------------------------------------------------

    int getNumCurveSegments() const
    {
        return m_CurveSegments.size();
    }

    const CurveSegment& getCurveSegment(CurveSegmentIndex ix) const
    {
        return m_CurveSegments[ix];
    }

    // Get body to which this cable's origin point is rigidly attached to.
    const Mobod& getOriginBody() const;

    // Get body to which this cable's termination point is rigidly attached to.
    const Mobod& getTerminationBody() const;

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

    Real getSurfaceConstraintTolerance() const
    {
        return m_ConstraintTolerance;
    }
    void setSurfaceConstraintTolerance(Real tolerance)
    {
        m_ConstraintTolerance = tolerance;
    }

    int getSurfaceProjectionMaxIter() const
    {
        return m_ProjectionMaxIter;
    }
    void setSurfaceProjectionMaxIter(int maxIter)
    {
        m_ProjectionMaxIter = maxIter;
    }

    Real getIntegratorAccuracy() const
    {
        return m_IntegratorAccuracy;
    }
    void setIntegratorAccuracy(Real accuracy)
    {
        m_IntegratorAccuracy = accuracy;
    }

    int getPathMaxIter() const
    {
        return m_PathMaxIter;
    }
    void setPathMaxIter(int maxIter)
    {
        m_PathMaxIter = maxIter;
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

//------------------------------------------------------------------------------

    // Find the number of CurveSegments that are in contact with an obstacle's
    // surface.
    size_t countActive(const State& s) const;

    int calcPathPoints(
        const State& state,
        std::function<void(Vec3 point_G)>& sink,
        int nPointsPerCurveSegment) const;

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
        return *m_Subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        return *m_Subsystem;
    }

    PosInfo& updPosInfo(const State& s) const;
    VelInfo& updVelInfo(const State& state) const;

    void calcPosInfo(const State& s, PosInfo& posInfo) const;
    void calcVelInfo(const State& s, VelInfo& velInfo) const;

//------------------------------------------------------------------------------
//                      HELPER FUNCTIONS
//------------------------------------------------------------------------------

    // Find the last contact point before the given curve segment, skipping
    // over any that are not in contact with their respective obstacle's
    // surface.
    Vec3 findPrevPoint(const State& state, CurveSegmentIndex ix) const;
    Vec3 findPrevPoint(const State& state, const CurveSegment& curve) const
    {
        return findPrevPoint(state, curve.getImpl().getIndex());
    }
    // Similarly find the first contact point after the given curve segment.
    Vec3 findNextPoint(const State& state, CurveSegmentIndex ix) const;
    Vec3 findNextPoint(const State& state, const CurveSegment& curve) const
    {
        return findNextPoint(state, curve.getImpl().getIndex());
    }
    // Similarly find the first curve segment before the given curve segment.
    const CurveSegment* findPrevActiveCurveSegment(
        const State& s,
        CurveSegmentIndex ix) const;
    // Similarly find the first curve segment after the given curve segment.
    const CurveSegment* findNextActiveCurveSegment(
        const State& s,
        CurveSegmentIndex ix) const;

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

    // Compute the straight line segments of this cable.
    Real calcLineSegments(
        const State& s,
        Vec3 p_O,
        Vec3 p_I,
        std::vector<LineSegment>& lines) const;

    // Call provided function for each curveSegment that is in contact with
    // their respective obstacle's surface.
    void callForEachActiveCurveSegment(
        const State& s,
        std::function<void(const CurveSegment::Impl&)> f) const;

//------------------------------------------------------------------------------

    // Reference back to the subsystem.
    CableSubsystem* m_Subsystem; // TODO just a pointer?

    MobilizedBodyIndex m_OriginBody;
    Vec3 m_OriginPoint;

    MobilizedBodyIndex m_TerminationBody;
    Vec3 m_TerminationPoint;

    Array_<CurveSegment, CurveSegmentIndex> m_CurveSegments{};

    Real m_PathAccuracy  = 1e-4;
    size_t m_PathMaxIter = 50; // TODO set to something reasonable.

    // For each curve segment the max allowed radial curvature.
    Real m_MaxCorrectionStepDeg = 10.; // TODO describe

    size_t m_ProjectionMaxIter = 50; // TODO set to reasonable value
    Real m_ConstraintTolerance = 1e-6;

    Real m_IntegratorAccuracy = 1e-6;

    // TOPOLOGY CACHE (set during realizeTopology())
    CacheEntryIndex m_PosInfoIx;
    CacheEntryIndex m_VelInfoIx;

    friend CurveSegment::Impl;
    friend CableSubsystemTestHelper;
};

} // namespace SimTK

#endif
