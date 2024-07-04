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

#include "Simbody.h"
#include "simbody/internal/CableSpan.h"
#include <iostream>

using namespace SimTK;

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
            const Real opacity = isInContactWithSurface ? 1. : 0.25;
            Transform X_GB =
                m_mbs->getMatterSubsystem()
                    .getMobilizedBody(m_cable.getObstacleMobilizedBodyIndex(ix))
                    .getBodyTransform(state);
            const Transform X_GS =
                X_GB.compose(m_cable.getObstacleXformSurfaceToBody(ix));

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

            if (!isInContactWithSurface) {
                continue;
            }

            // Draw the frenet frames at the geodesic boundary points.
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

    // Mobilizer for path origin.
    MobilizedBody::Translation cableOriginBody(
        matter.Ground(),
        Vec3(-8., 0.1, 0.),
        aBody,
        Transform());

    // Mobilizer for path termination.
    MobilizedBody::Translation cableTerminationBody(
        matter.Ground(),
        Transform(Vec3(10., 1.0, -1.)),
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

    // Add torus obstacle.
    cable.addObstacle(
        matter.Ground(), // Obstacle mobilizer body.
        Transform(
            Rotation(0.5 * Pi, YAxis),
            Vec3{-4., 0., 0.}), // Surface to body transform.
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Torus(1., 0.2)), // Obstacle geometry.
        {0.1, 0.2, 0.}                            // Initial contact point hint.
    );

    // Add ellipsoid obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Vec3{-2., 0., 0.}),
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Ellipsoid({1.5, 2.6, 1.})),
        {0.0, 0., 1.1});

    // Add sphere obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Vec3{2., 0., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(1.)),
        {0.1, 1.1, 0.});

    const bool show = false;
    system.setUseUniformBackground(true); // no ground plane in display
    // Visualize the system.
    std::unique_ptr<Visualizer> viz(show? new Visualizer(system) : nullptr);

    if (viz) {
        viz->setShowFrameNumber(true);
        viz->addDecorationGenerator(new CableDecorator(system, cable));
    }

    // Initialize the system and s.
    system.realizeTopology();
    State s = system.getDefaultState();

    Real angle = 0.;
    for (int i = 0; i < 150; ++i) {
        system.realize(s, Stage::Report);

        for (CableSpanObstacleIndex ix(0); ix < cable.getNumObstacles(); ++ix) {
            SimTK_ASSERT1_ALWAYS(
                cable.isInContactWithObstacle(s, ix),
                "lost contact with obstacle %d",
                ix);
        }

        if(viz) {
            viz->report(s);
        }
        cable.storeCurrentPath(s);

        // Move the cable origin.
        angle += 0.01;
        cableOriginBody.setQ(
            s,
            Vec3(sin(angle), 5. * sin(angle * 1.5), 5. * sin(angle * 2.)));
    }
}
