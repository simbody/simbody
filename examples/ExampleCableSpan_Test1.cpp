/*-----------------------------------------------------------------------------
                Simbody(tm) Example: Cable Over Smooth Surfaces
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

/* This example is for experimenting with a CableSpan over obstacle surfaces.
The cable endpoints are manually repositioned to control the next path
solving problem. */

#include "SimTKcommon/internal/Timing.h"
#include "Simbody.h"
#include "simbody/internal/CableSpan.h"
#include <cstdio>
#include <iostream>

using namespace SimTK;

// A helper class for printing interesting cable outputs.
class ShowStuff : public PeriodicEventReporter {
public:
    ShowStuff(const CableSubsystem& cables, Real interval) :
        PeriodicEventReporter(interval), m_cables(cables)
    {}

    static void showHeading(std::ostream& o)
    {
        printf(
            "%10s %10s %10s %10s",
            "Cable index",
            "iterations",
            "path error",
            "length");
    }

    /** This is the implementation of the EventReporter virtual. **/
    void handleEvent(const State& s) const override
    {
        for (CableSpanIndex cableIx(0); cableIx < m_cables.getNumCables();
             ++cableIx) {
            const CableSpan& cable = m_cables.getCable(cableIx);
            printf(
                "%d %d %g %g",
                static_cast<int>(cableIx),
                cable.getNumSolverIterations(s),
                cable.getSmoothness(s),
                cable.calcLength(s));
        }
    }

private:
    const CableSubsystem& m_cables;
};

/* A helper class for drawing interesting things of a CableSpan. */
class CableDecorator : public SimTK::DecorationGenerator {
public:
    CableDecorator(MultibodySystem& mbs, const CableSpan& cable) :
        m_mbs(&mbs), m_cable(cable)
    {
        for (CableSpanObstacleIndex ix(0); ix < m_cable.getNumObstacles();
             ++ix) {
            m_obstacleDecorations.push_back(
                m_cable.getObstacleContactGeometry(ix)
                    .createDecorativeGeometry()
                    .setResolution(3));
            m_obstacleDecorationsOffsets.push_back(
                m_obstacleDecorations.back().getTransform());
        }
    }

    virtual void generateDecorations(
        const State& state,
        Array_<DecorativeGeometry>& decorations) override
    {
        for (CableSpanObstacleIndex ix(0); ix < m_cable.getNumObstacles();
             ++ix) {
            // Draw the obstacle surface.
            // Green if cable is in contact with surface, grey otheriwse.
            const ContactGeometry& geometry =
                m_cable.getObstacleContactGeometry(ix);
            const bool isInContactWithSurface =
                m_cable.isInContactWithObstacle(state, ix);
            const Vec3 color   = isInContactWithSurface ? Yellow : Gray;
            const Real opacity = isInContactWithSurface ? 0.5 : 0.25;
            Transform X_GB =
                m_mbs->getMatterSubsystem()
                    .getMobilizedBody(m_cable.getObstacleMobilizedBodyIndex(ix))
                    .getBodyTransform(state);
            const Transform X_GS =
                X_GB.compose(m_cable.getObstacleTransformSurfaceToBody(ix));

            const Transform X_SD = m_obstacleDecorationsOffsets.at(ix);
            const Transform X_GD = X_GS.compose(X_SD);
            // Draw the obstacle's local frame.
            // This is the frame that you define the contact point hint in.
            decorations.push_back(
                DecorativeFrame(0.5).setTransform(X_GS).setColor(Purple));
            // Draw the obstacle contact geometry.
            decorations.push_back(m_obstacleDecorations.at(ix)
                                      .setTransform(X_GD)
                                      .setColor(color)
                                      .setOpacity(opacity));

            // Draw the initial contact point hints (these are user-defined) as
            // an orange line with a point.
            const Vec3 x_PS = m_cable.getObstacleContactPointHint(ix);
            decorations.push_back(
                DecorativeLine(X_GS.p(), X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Orange)
                    .setLineThickness(3));
            decorations.push_back(
                DecorativePoint(X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Orange));

            // Only draw curve segment related decorations if the cable is in
            // contact with the corresponding obstacle.
            if (!isInContactWithSurface) {
                continue;
            }

            // Draw a fixed number of sphere points per curve segment.
            /* const int numDecorativePoints = 4; */
            /* m_cable.calcCurveSegmentResampledPoints( */
            /*     state, */
            /*     ix, */
            /*     numDecorativePoints, */
            /*     [&](Vec3 x_G) */
            /*     { */
            /*         decorations.push_back( */
            /*             DecorativeSphere(0.1).setTransform(x_G).setColor(Blue)); */
            /*     }); */

            // Draw the Frenet frames at the geodesic boundary points.
            decorations.push_back(DecorativeFrame(0.2).setTransform(
                m_cable.calcCurveSegmentInitialFrenetFrame(state, ix)));
            decorations.push_back(DecorativeFrame(0.2).setTransform(
                m_cable.calcCurveSegmentFinalFrenetFrame(state, ix)));
        }
    }

    MultibodySystem* m_mbs;
    CableSpan m_cable;
    int m_lineThickness = 3;

    Array_<DecorativeGeometry, CableSpanObstacleIndex> m_obstacleDecorations;
    Array_<Transform, CableSpanObstacleIndex> m_obstacleDecorationsOffsets;
};

int main()
{
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // A dummy body.
    Body::Rigid aBody(MassProperties(1.0, Vec3(0), Inertia(1)));

    // Offset between origin & termination.
    const Real offset = 0.1;

    // Mobilizer for path origin.
    MobilizedBody::Translation cableOriginBody(
        matter.Ground(),
        Transform(Vec3(-offset, 0., 0.)),
        aBody,
        Transform());

    // Mobilizer for path termination.
    MobilizedBody::Translation cableTerminationBody(
        matter.Ground(),
        Transform(Vec3(offset, 0., 0.)),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        cableOriginBody,
        Vec3{0.},
        cableTerminationBody,
        Vec3{0.});

    // Add some obstacles to the cable.
    cable.addObstacle(
        matter.Ground(),
        Transform(
            Rotation(0.5 * Pi, YAxis),
            Vec3{0., 4., 0.}),
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Torus(2., 0.5)),
        {0., -0.1, 0.});

    cable.setCurveSegmentAccuracy(1e-10);

    // Visaulize the system.
    system.setUseUniformBackground(true); // no ground plane in display
    Visualizer viz(system);
    viz.addDecorationGenerator(new CableDecorator(system, cable));
    ShowStuff showStuff(cables, 1e-3);

    // Initialize the system and s.
    system.realizeTopology();
    State s = system.getDefaultState();

    system.realize(s, Stage::Report);
    viz.report(s);
    cable.calcLength(s);
    showStuff.handleEvent(s);
    std::cout << "solved in " << cable.getNumSolverIterations(s) << " iterations\n";
}
