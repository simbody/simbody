#ifndef SimTK_SIMBODY_CABLE_SPAN_SUBSYSTEM_IMPL_H_
#define SimTK_SIMBODY_CABLE_SPAN_SUBSYSTEM_IMPL_H_

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

#include "simbody/internal/CableSpan.h"
#include "simbody/internal/MultibodySystem.h"
#include "CableSpan_Impl.h"
#include "CableSpan_CurveSegment_Impl.h"

namespace SimTK
{

//==============================================================================
//                         SUBSYSTEM :: IMPL
//==============================================================================
class CableSubsystem::Impl : public Subsystem::Guts
{
public:
    /** This is a helper struct that is used by a CableSpan to compute the
    position level cache entry.
    After computing the position level cache, this data is discarded. */
    struct SolverData
    {
        // Number of degrees of freedom of a geodesic.
        static constexpr int GEODESIC_DOF = 4;

        // Number of path error constraints per curve segment.
        static constexpr int NUMBER_OF_CONSTRAINTS = 4;

        // The matrix dimensions and number of straight line segments are
        // determined by the number of CurveSegments that are in contact with
        // their respective obstacle's surface.
        SolverData(int nActive) : nCurves(nActive)
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
        // TODO Best solver?
        FactorQTZ inverse; // TODO rename to inv
        int nCurves = -1;
    };

    // Cache entry for holding the SolverData of different dimensions.
    struct CacheEntry
    {
        CacheEntry() = default;
        SolverData& updOrInsert(int nActive)
        {
            if (nActive <= 0) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Cannot produce solver data of zero dimension");
            }

            for (int i = m_Data.size(); i < nActive; ++i) {
                m_Data.emplace_back(i + 1);
            }

            return m_Data.at(nActive - 1);
        }

        // Solver data of suitable for solving the cable path given a number of
        // active CurveSegments. The ith element in the vector is used for
        // solving a path with i+1 active CurveSegments.
        std::vector<SolverData> m_Data;
    };

    // TODO fix rule of five
    Impl()
    {}
    ~Impl()
    {}

    void realizeTopology(State& state)
    {
        CacheEntry cache{};
        m_CacheIx = allocateCacheEntry(
            state,
            Stage::Instance,
            Stage::Infinity,
            new Value<CacheEntry>(cache));
    }

    CableSpanIndex adoptCable(
            CableSubsystem& subsystemHandle,
            CableSpan& cable)
    {
        invalidateSubsystemTopologyCache();

        CableSpanIndex cableIx(cables.size());
        cable.updImpl().setSubsystem(subsystemHandle, cableIx);
        cables.push_back(cable);
        return cableIx;
    }

    int getNumCables() const
    {
        return cables.size();
    }

    const CableSpan& getCable(CableSpanIndex index) const
    {
        return cables[index];
    }

    CableSpan& updCable(CableSpanIndex index)
    {
        return cables[index];
    }

    // Return the MultibodySystem which owns this WrappingPathSubsystem.
    const MultibodySystem& getMultibodySystem() const
    {
        return MultibodySystem::downcast(getSystem());
    }

    // Return the SimbodyMatterSubsystem from which this WrappingPathSubsystem
    // gets the bodies to track.
    const SimbodyMatterSubsystem& getMatterSubsystem() const
    {
        return getMultibodySystem().getMatterSubsystem();
    }

    CacheEntry& updSolverData(const State& state) const
    {
        return Value<CacheEntry>::updDowncast(updCacheEntry(state, m_CacheIx));
    }

    SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
    Impl* cloneImpl() const override
    {
        return new Impl(*this);
    }

    int calcDecorativeGeometryAndAppendImpl(
        const State& state,
        Stage stage,
        Array_<DecorativeGeometry>& decorations) const override;

    // Allocate state variables.
    int realizeSubsystemTopologyImpl(State& state) const override;

    // TOPOLOGY STATE
    Array_<CableSpan, CableSpanIndex> cables;

    CacheEntryIndex m_CacheIx;

    friend CableSubsystemTestHelper;
};

} // namespace SimTK

#endif
