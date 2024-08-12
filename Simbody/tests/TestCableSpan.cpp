/*-----------------------------------------------------------------------------
                Simbody(tm) Test: Cable Over Several Smooth Surfaces
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

/**
This file contains tests for the CableSpan's path over different obstacles.
**/

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

    void generateDecorations(
        const State& state,
        Array_<DecorativeGeometry>& decorations) override
    {
        for (CableSpanObstacleIndex ix(0); ix < m_cable.getNumObstacles();
             ++ix) {
            // Draw the obstacle surface.
            const ContactGeometry& geometry =
                m_cable.getObstacleContactGeometry(ix);
            // If cable is not in contact with the surface grey it out.
            const bool isInContactWithSurface =
                m_cable.isInContactWithObstacle(state, ix);
            const Vec3 color   = isInContactWithSurface ? Yellow : Gray;
            const Real opacity = isInContactWithSurface ? 0.5 : 0.25;
            // Transform from Ground to obstacle body.
            Transform X_GB =
                m_mbs->getMatterSubsystem()
                    .getMobilizedBody(m_cable.getObstacleMobilizedBodyIndex(ix))
                    .getBodyTransform(state);
            // Transform from Ground to obstacle contact surface offset frame.
            const Transform X_GS =
                X_GB.compose(m_cable.getObstacleXformSurfaceToBody(ix));
            // Transform from ground to decoration surface.
            const Transform X_GD =
                X_GS.compose(m_obstacleDecorationsOffsets.at(ix));
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
            // a line and a point.
            const Vec3 x_PS = m_cable.getObstacleContactPointHint(ix);
            decorations.push_back(
                DecorativeLine(X_GS.p(), X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Green)
                    .setLineThickness(3));
            decorations.push_back(
                DecorativePoint(X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Green));
        }
    }

    MultibodySystem* m_mbs;
    CableSpan m_cable;
    Array_<DecorativeGeometry, CableSpanObstacleIndex> m_obstacleDecorations;
    Array_<Transform, CableSpanObstacleIndex> m_obstacleDecorationsOffsets;
};

/** Simple CableSpon path with known solution.

Wrap a cable over (in order):
1. Torus
2. Ellipsoid
3. Sphere
4. Cylinder

We wrap the cable conveniently over the obstacles such that each curve
segment becomes a circular-arc shape. This allows us to check the results by
hand. **/
void testSimpleCable()
{
    const bool show = false;

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // A dummy body.
    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    // Mobilizer for path origin.
    MobilizedBody::Translation cableOriginBody(
        matter.Ground(),
        Vec3(0.),
        aBody,
        Transform());

    // Mobilizer for path termination.
    MobilizedBody::Translation cableTerminationBody(
        matter.Ground(),
        Transform(Vec3(-4., 0., 0.)),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        cableOriginBody,
        Vec3{0.},
        cableTerminationBody,
        Vec3{0.});
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-6);

    // Add initial torus obstacle.
    MobilizedBody::Translation torusBody(
        matter.Ground(),
        Transform(Vec3(0., 10., 0.)),
        aBody,
        Transform());
    cable.addObstacle(
        torusBody,
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{1., 1., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Torus(10., 1.)),
        {0., -9., 0.});

    // Add ellipsoid obstacle.
    MobilizedBody::Translation ellipsoidBody(
        torusBody,
        Transform(Vec3(4., -10., 0.)),
        aBody,
        Transform());
    cable.addObstacle(
        ellipsoidBody,
        Transform(Vec3{0.5, 1., 0.}),
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Ellipsoid({1., 1., 0.75})),
        {1., 1., 0.});

    // Add sphere obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Vec3{4., -2., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(1.5)),
        {0., -1., 0.});

    // Add cylinder obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Vec3{-2., -1.5, 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Cylinder(2.)),
        Vec3{0., -1., 0.});

    // Optionally visualize the system.
    system.setUseUniformBackground(true); // no ground plane in display
    std::unique_ptr<Visualizer> viz(show ? new Visualizer(system) : nullptr);

    if (viz) {
        viz->setShowFrameNumber(true);
        viz->addDecorationGenerator(new CableDecorator(system, cable));
    }

    // Initialize the system and state.
    system.realizeTopology();
    State s = system.getDefaultState();

    // Compute the CableSpan's path.
    system.realize(s, Stage::Report);
    const Real cableLength = cable.calcLength(s);
    if (viz) {
        viz->report(s);
    }

    SimTK_ASSERT_ALWAYS(
        cable.getNumObstacles() == 4,
        "Invalid number of obstacles");

    SimTK_ASSERT2_ALWAYS(
        cable.getSmoothness(s) <= cable.getSmoothnessTolerance(),
        "Test failed: Cable smoothness (=%e) must be smaller than set tolerance (=%e)",
        cable.getSmoothness(s),
        cable.getSmoothnessTolerance());

    std::array<std::string, 4> obsNames{
        "torus",
        "ellipsoid",
        "sphere",
        "cylinder"};

    // Note that the length deviates because the path is solved up to angle
    // tolerance, not length tolerance.
    const Real lengthTolerance = 1e-5;
    auto assertCurveSegmentLength =
        [&](CableSpanObstacleIndex obsIx, Real obsRadius)
    {
        SimTK_ASSERT1_ALWAYS(
            cable.isInContactWithObstacle(s, obsIx),
            "Cable not in contact with %s obstacle",
            obsNames.at(obsIx).c_str());

        const Real angle          = 0.5 * Pi;
        const Real expectedLength = angle * obsRadius;
        const Real gotLength      = cable.calcCurveSegmentArcLength(s, obsIx);

        SimTK_ASSERT4_ALWAYS(
            std::abs(gotLength - expectedLength) < lengthTolerance,
            "%s curve segment length (=%f) does not match expected length (=%f), with error %e",
            obsNames.at(obsIx).c_str(),
            gotLength,
            expectedLength,
            gotLength - expectedLength);
    };

    assertCurveSegmentLength(CableSpanObstacleIndex(0), 1.);
    assertCurveSegmentLength(CableSpanObstacleIndex(1), 1.);
    assertCurveSegmentLength(CableSpanObstacleIndex(2), 1.5);
    assertCurveSegmentLength(CableSpanObstacleIndex(3), 2.);

    // Sum all straight line segment lengths + curve line segment lengths.
    const Real sumStraightLineSegmentLengths = 1. + 3.5 + 3. + 6. + 1.5;
    const Real sumCurveLineSegmentLengths    = 0.5 * Pi * (1. + 1. + 1.5 + 2.);
    const Real expectedTotalCableLength =
        sumStraightLineSegmentLengths + sumCurveLineSegmentLengths;
    const Real gotTotalCableLength = cable.calcLength(s);
    SimTK_ASSERT3_ALWAYS(
        std::abs(expectedTotalCableLength - gotTotalCableLength) <
            lengthTolerance,
        "Expected cable length (=%f) does not match computed cable length (=%f), error = %e",
        expectedTotalCableLength,
        gotTotalCableLength,
        expectedTotalCableLength - gotTotalCableLength);

    // We should pass the generic tests:
    CableSubsystemTestHelper().testCurrentPath(s, cables, std::cout);

    const Real frenetFrameTolerance     = 1e-6;
    auto assertCurveSegmentFrenetFrames = [&](CableSpanObstacleIndex obsIx,
                                              const Transform& expected_X_GP,
                                              const Transform& expected_X_GQ)
    {
        const Transform& got_X_GP =
            cable.calcCurveSegmentInitialFrenetFrame(s, obsIx);
        const Transform& got_X_GQ =
            cable.calcCurveSegmentFinalFrenetFrame(s, obsIx);

        SimTK_ASSERT1_ALWAYS(
            (expected_X_GP.p() - got_X_GP.p()).norm() < frenetFrameTolerance,
            "%s curve segment position at initial contact point incorrect",
            obsNames.at(obsIx).c_str());
        SimTK_ASSERT1_ALWAYS(
            (expected_X_GP.R().asMat33() - got_X_GP.R().asMat33()).norm() <
                frenetFrameTolerance,
            "%s curve segment frame orientation at initial contact point incorrect",
            obsNames.at(obsIx).c_str());

        SimTK_ASSERT1_ALWAYS(
            (expected_X_GQ.p() - got_X_GQ.p()).norm() < frenetFrameTolerance,
            "%s curve segment position at final contact point incorrect",
            obsNames.at(obsIx).c_str());
        SimTK_ASSERT1_ALWAYS(
            (expected_X_GQ.R().asMat33() - got_X_GQ.R().asMat33()).norm() <
                frenetFrameTolerance,
            "%s curve segment frame orientation at final contact point incorrect",
            obsNames.at(obsIx).c_str());
    };

    assertCurveSegmentFrenetFrames(
        CableSpanObstacleIndex(0),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(0., 1., 0.)),
                XAxis,
                UnitVec3(Vec3(-1., 0., 0.)),
                YAxis),
            Vec3(0., 1., 0.)),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(1., 0., 0.)),
                XAxis,
                UnitVec3(Vec3(0., 1., 0.)),
                YAxis),
            Vec3(1., 2., 0.)));
    assertCurveSegmentFrenetFrames(
        CableSpanObstacleIndex(1),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(1., 0., 0.)),
                XAxis,
                UnitVec3(Vec3(0., 1., 0.)),
                YAxis),
            Vec3(4.5, 2., 0.)),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(0., -1., 0.)),
                XAxis,
                UnitVec3(Vec3(1., 0., 0.)),
                YAxis),
            Vec3(5.5, 1., 0.)));
    assertCurveSegmentFrenetFrames(
        CableSpanObstacleIndex(2),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(0., -1., 0.)),
                XAxis,
                UnitVec3(Vec3(1., 0., 0.)),
                YAxis),
            Vec3(5.5, -2., 0.)),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(-1., 0., 0.)),
                XAxis,
                UnitVec3(Vec3(0., -1., 0.)),
                YAxis),
            Vec3(4., -3.5, 0.)));
    assertCurveSegmentFrenetFrames(
        CableSpanObstacleIndex(3),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(-1., 0., 0.)),
                XAxis,
                UnitVec3(Vec3(0., -1., 0.)),
                YAxis),
            Vec3(-2., -3.5, 0.)),
        Transform(
            Rotation().setRotationFromTwoAxes(
                UnitVec3(Vec3(0., 1., 0.)),
                XAxis,
                UnitVec3(Vec3(-1., 0., 0.)),
                YAxis),
            Vec3(-4., -1.5, 0.)));

    Vector_<SpatialVec> expectedBodyForcesInG(5, {{0., 0., 0.}, {0., 0., 0.}});
    // Ground body - Sphere and cylinder are connected to it.
    expectedBodyForcesInG[0] = (SpatialVec{{0., 0., 1.5}, {0., 2., 0.}});
    // Cable origin body.
    expectedBodyForcesInG[1] = (SpatialVec{{0., 0., 0.}, {0., -1., 0.}});
    // Cable termination body.
    expectedBodyForcesInG[2] = (SpatialVec{{0., 0., 0.}, {0., -1., 0.}});
    // Torus body.
    expectedBodyForcesInG[3] = (SpatialVec{{0., 0., 8.}, {1., -1., 0.}});
    // Ellipsoid body.
    expectedBodyForcesInG[4] = (SpatialVec{{0., 0., 0.5}, {-1., -1., 0.}});

    const Real tension = Pi;
    Vector_<SpatialVec> gotBodyForcesInG(5, {{0., 0., 0.}, {0., 0., 0.}});
    cable.applyBodyForces(s, tension, gotBodyForcesInG);

    for (int i = 0; i < gotBodyForcesInG.size(); ++i) {
        const Real tolerance = 1e-5;
        SimTK_ASSERT1_ALWAYS(
            (expectedBodyForcesInG[i] * tension - gotBodyForcesInG[i]).norm() <
                tolerance,
            "Unexpected force applied to body %i",
            i);
    }
}

/** Test computed cable path over all supported surfaces, testing geodesics,
jacobians, and kinematics.

The cable wraps over (in order):
1. Torus,
2. Ellipsoid,
3. Sphere,
4. Cylinder,
5. Bicubic patch,
6. Torus.

The flag assertCableLengthDerivative is used to switch between verifying that
the computed cable length derivative matches the change in length during
simulation, OR verifying that all computed geodesics and jacobians are
correct. Doing both requires too much time. **/
void testAllSurfaceKinds(bool assertCableLengthDerivative)
{
    const bool show = false;

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // A dummy body.
    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    // Mobilizer for path origin.
    MobilizedBody::Translation cableOriginBody(
        matter.Ground(),
        Vec3(-8., 0.1, 0.),
        aBody,
        Transform());

    // Mobilizer for path termination.
    MobilizedBody::Translation cableTerminationBody(
        matter.Ground(),
        Transform(Vec3(20., 1.0, -1.)),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        cableOriginBody,
        Vec3{0.},
        cableTerminationBody,
        Vec3{0.});
    cable.setCurveSegmentAccuracy(1e-9);
    cable.setSmoothnessTolerance(1e-4);

    // Add initial torus obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{-4., 0., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Torus(1., 0.2)),
        {0.1, 0.2, 0.});

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

    // Add cylinder obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Rotation(0.5 * Pi, XAxis), Vec3{5., 0., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Cylinder(1.)),
        Vec3{0., -1., 0.});

    // Add bicubic surface obstacle.
    {
        Real patchScaleX = 2.0;
        Real patchScaleY = 2.0;
        Real patchScaleF = 0.75;

        constexpr int Nx = 4;
        constexpr int Ny = 4;

        const Real xData[Nx] = {-2, -1, 1, 2};
        const Real yData[Ny] = {-2, -1, 1, 2};

        const Real fData[Nx * Ny] =
            {2, 3, 3, 1, 0, 1.5, 1.5, 0, 0, 1.5, 1.5, 0, 2, 3, 3, 1};

        const Vector x_(Nx, xData);
        const Vector y_(Ny, yData);
        const Matrix f_(Nx, Ny, fData);

        Vector x = patchScaleX * x_;
        Vector y = patchScaleY * y_;
        Matrix f = patchScaleF * f_;

        BicubicSurface patch(x, y, f, 0);
        Transform patchTransform(
            Rotation(0.5 * Pi, Vec3(0., 0., 1.)),
            Vec3(10., 0., -1.));

        cable.addObstacle(
            matter.Ground(),
            patchTransform,
            std::shared_ptr<const ContactGeometry>(
                new ContactGeometry::SmoothHeightMap(patch)),
            Vec3{0., 0., 1.});
    }

    // Add final torus obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{14., 0., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Torus(1., 0.2)),
        {0.1, 0.2, 0.});

    // Optionally visualize the system.
    system.setUseUniformBackground(true); // no ground plane in display
    std::unique_ptr<Visualizer> viz(show ? new Visualizer(system) : nullptr);

    if (viz) {
        viz->setShowFrameNumber(true);
        viz->addDecorationGenerator(new CableDecorator(system, cable));
    }

    // Initialize the system and state.
    system.realizeTopology();
    State s = system.getDefaultState();

    // Use this to assert cable length time derivative.
    Real prevCableLength = NaN;

    // Let the cable end points be parameterized by an angle, and draw the path
    // for different angles. If we want to check the CableSpan::lengthDot we
    // will use a smaller stepsize, and smaller final angle, and not run
    // CableSubsystemTestHelper on each computed path.
    const Real dAngle     = assertCableLengthDerivative ? 1e-4 : 0.05;
    const Real finalAngle = assertCableLengthDerivative ? 0.5 * Pi : 2. * Pi;
    for (Real angle = 0.; angle < finalAngle; angle += dAngle) {

        // Move the cable end points.
        cableOriginBody.setQ(
            s,
            Vec3(
                1.1 * sin(angle),
                5. * sin(angle * 1.5),
                5. * sin(angle * 2.)));
        cableOriginBody.setU(
            s,
            Vec3(
                1.1 * cos(angle),
                5. * 1.5 * cos(angle * 1.5),
                5. * 2. * cos(angle * 2.)));
        cableTerminationBody.setQ(
            s,
            Vec3(
                0.1 * sin(angle),
                2. * sin(angle * 0.7),
                5. * sin(angle * 1.3)));
        cableTerminationBody.setU(
            s,
            Vec3(
                0.1 * cos(angle),
                2. * 0.7 * cos(angle * 0.7),
                5. * 1.3 * cos(angle * 1.3)));

        // Compute the CableSpan's path.
        system.realize(s, Stage::Report);
        const Real cableLength = cable.calcLength(s);
        if (viz) {
            viz->report(s);
        }

        // Check that the geodesics and path error vector & jacobian are
        // correct.
        if (!assertCableLengthDerivative) {
            CableSubsystemTestHelper().testCurrentPath(s, cables, std::cout);
        }

        // Assert length derivative using the change in length.
        if (assertCableLengthDerivative && !isNaN(prevCableLength)) {
            const Real tolerance = 5e-3;
            const Real expectedCableLengthDot =
                (cableLength - prevCableLength) / dAngle;
            const Real gotCableLengthDot = cable.calcLengthDot(s);
            SimTK_ASSERT4_ALWAYS(
                std::abs(gotCableLengthDot - expectedCableLengthDot) <
                    tolerance,
                "Test failed: Cable length dot (=%f) does not match expected cable length dot (=%f), with error %f, at angle %f",
                gotCableLengthDot,
                expectedCableLengthDot,
                gotCableLengthDot - expectedCableLengthDot,
                angle);
        }

        // Total cable length should be longer than direct distance between the
        // cable endpoints.
        const Real distanceBetweenEndPoints =
            (cableTerminationBody.getBodyOriginLocation(s) -
             cableOriginBody.getBodyOriginLocation(s))
                .norm();
        SimTK_ASSERT2_ALWAYS(
            cableLength > distanceBetweenEndPoints,
            "Test failed: Cable length (=%f) smaller than distance between end points (=%f)",
            cableLength,
            distanceBetweenEndPoints);

        // Make sure that we acutally solved the path up to tolerance.
        SimTK_ASSERT2_ALWAYS(
            cable.getSmoothness(s) <= cable.getSmoothnessTolerance(),
            "Test failed: Cable smoothness (=%e) must be smaller than set tolerance (=%e)",
            cable.getSmoothness(s),
            cable.getSmoothnessTolerance());

        prevCableLength = cableLength;
    }
    std::cout
        << "PASSED TEST: testAllSurfaceKinds (assertCableLengthDerivative = "
        << assertCableLengthDerivative << std::endl;
}

/** CableSpan touchdown and liftoff test.

Create a simple touchdown and liftoff case on a:
- torus,
- ellipsoid,
- sphere,
- cylinder.

A cable is spanned over each obstacle individually, so there are 4 cables, each
with one obstacle. **/
void testTouchdownAndLiftoff()
{
    const bool show = false;

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // A dummy body.
    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    // Mobilizer for path origin.
    MobilizedBody::Translation cableOriginBody(
        matter.Ground(),
        Vec3(-2., 0., 0.),
        aBody,
        Transform());

    // Mobilizer for path termination.
    MobilizedBody::Translation cableTerminationBody(
        matter.Ground(),
        Transform(Vec3(2., 0., 0.)),
        aBody,
        Transform());

    auto createCable = [&]()
    {
        CableSpan cable(
            cables,
            cableOriginBody,
            Vec3{0.},
            cableTerminationBody,
            Vec3{0.});
        cable.setCurveSegmentAccuracy(1e-12);
        cable.setSmoothnessTolerance(1e-6);
        return cable;
    };

    // Create a cable with a torus obstacle.
    {
        CableSpan cable = createCable();
        cable.addObstacle(
            matter.Ground(),
            Transform(Rotation(0.5 * Pi, YAxis), Vec3{0., 1.75, 0.}),
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Torus(2., 0.25)));
    }

    // Create a cable with an ellipsoid obstacle.
    {
        CableSpan cable = createCable();
        cable.addObstacle(
            matter.Ground(),
            Transform(Vec3{0., -0.5, 0.}),
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Ellipsoid({1., 0.5, 0.75})));
    }

    // Create a cable with a sphere obstacle.
    {
        CableSpan cable = createCable();
        cable.addObstacle(
            matter.Ground(),
            Transform(Vec3{0., -1.5, 0.}),
            std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(1.5)));
    }

    // Create a cable with a cylinder obstacle.
    {
        CableSpan cable = createCable();
        // Add cylinder obstacle.
        cable.addObstacle(
            matter.Ground(),
            Transform(Rotation(0. * Pi, XAxis), Vec3{0., -2., 0.}),
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Cylinder(2.)));
    }

    // Optionally visualize the system.
    system.setUseUniformBackground(true); // no ground plane in display
    std::unique_ptr<Visualizer> viz(show ? new Visualizer(system) : nullptr);
    if (viz) {
        viz->setShowFrameNumber(true);
        for (CableSpanIndex cableIx(0); cableIx < cables.getNumCables();
             ++cableIx) {
            viz->addDecorationGenerator(
                new CableDecorator(system, cables.getCable(cableIx)));
        }
    }

    // Initialize the system and state.
    system.realizeTopology();
    State s = system.getDefaultState();

    SimTK_ASSERT1_ALWAYS(
        cables.getNumCables() == 4,
        "Unexpected number of cables (=%i)",
        cables.getNumCables());

    for (Real angle = 1e-2; angle < 4. * Pi; angle += 0.02) {
        const Real yCoord = 0.1 * sin(angle);

        // Move the cable end points.
        cableOriginBody.setQ(s, Vec3(0., yCoord, 0.));
        cableTerminationBody.setQ(s, Vec3(0., yCoord, 0.));

        system.realize(s, Stage::Report);

        // All obstacles are positioned such that for negative yCoord of the
        // endpoints, the cable touches down on the obstacle.
        for (CableSpanIndex cableIx(0); cableIx < cables.getNumCables();
             ++cableIx) {
            cables.getCable(cableIx).calcLength(s);
            const bool gotContactStatus =
                cables.getCable(cableIx).isInContactWithObstacle(
                    s,
                    CableSpanObstacleIndex(0));
            const bool expectedContactStatus = yCoord < 0.;
            SimTK_ASSERT4_ALWAYS(
                gotContactStatus == expectedContactStatus,
                "Cable %i expected contact status (=%i) does not match computed contact status (=%i) at yCoord = %f",
                cableIx,
                expectedContactStatus,
                gotContactStatus,
                angle);
        }

        if (viz) {
            viz->report(s);
        }
    }
}

int main()
{
    testSimpleCable();
    testTouchdownAndLiftoff();
    testAllSurfaceKinds(false); // Test all geodesics and jacobians.
    testAllSurfaceKinds(true);  // Test length derivative.
}
