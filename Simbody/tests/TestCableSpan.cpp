/*-----------------------------------------------------------------------------
            Simbody(tm) Test: Cable Over Smooth Surfaces and Via Points
-------------------------------------------------------------------------------
 Copyright (c) 2024 Authors.
 Authors: Pepijn van den Bos
 Contributors: Nicholas Bianco

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

using namespace SimTK;

/**
This file contains tests for the CableSpan's path over different obstacles and
via points.
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
                X_GB.compose(m_cable.getObstacleTransformSurfaceToBody(ix));
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

        for (CableSpanViaPointIndex ix(0); ix < m_cable.getNumViaPoints();
                ++ix) {
            // Draw the via point.
            const Vec3 x_G = m_cable.calcViaPointLocation(state, ix);
            decorations.push_back(
                DecorativeSphere(0.1)
                    .setTransform(x_G)
                    .setRepresentation(DecorativeGeometry::DrawWireframe)
                    .setColor(Cyan)
                    .setOpacity(0.1));
        }
    }

    MultibodySystem* m_mbs;
    CableSpan m_cable;
    Array_<DecorativeGeometry, CableSpanObstacleIndex> m_obstacleDecorations;
    Array_<Transform, CableSpanObstacleIndex> m_obstacleDecorationsOffsets;
};

/** Simple CableSpan path with known solution.

Wrap a cable over (in order):
1. Torus
2. Ellipsoid
3. Torus
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
    cable.setSmoothnessTolerance(1e-7);

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
            new ContactGeometry::Ellipsoid({1., 1., 6.})),
        {1., 1., 0.5});

    // Add another torus obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{4., -2. - 10., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Torus(10., 1.5)),
        {0., 1., 0.});

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
    cable.calcLength(s);
    if (viz) {
        viz->report(s);
    }

    SimTK_ASSERT_ALWAYS(
        cable.getNumObstacles() == 4,
        "Invalid number of obstacles");

    SimTK_ASSERT_ALWAYS(
        cable.getNumViaPoints() == 0,
        "Invalid number of via points");

    SimTK_ASSERT2_ALWAYS(
        cable.getSmoothness(s) <= cable.getSmoothnessTolerance(),
        "Test failed: Cable smoothness (=%e) must be smaller than set "
        "tolerance (=%e)",
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
            "%s curve segment length (=%f) does not match expected length "
            "(=%f), with error %e",
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
        "Expected cable length (=%f) does not match computed cable length "
        "(=%f), error = %e",
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
            "%s curve segment frame orientation at initial contact point "
            "incorrect",
            obsNames.at(obsIx).c_str());

        SimTK_ASSERT1_ALWAYS(
            (expected_X_GQ.p() - got_X_GQ.p()).norm() < frenetFrameTolerance,
            "%s curve segment position at final contact point incorrect",
            obsNames.at(obsIx).c_str());
        SimTK_ASSERT1_ALWAYS(
            (expected_X_GQ.R().asMat33() - got_X_GQ.R().asMat33()).norm() <
                frenetFrameTolerance,
            "%s curve segment frame orientation at final contact point "
            "incorrect",
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
    expectedBodyForcesInG[1] = (SpatialVec{{0., 0., 0.}, {0., 1., 0.}});
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

/** Simple CableSpan comprised of via points with known solution. **/
void testViaPoints()
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

    // First via point.
    MobilizedBody::Translation cableViaPoint1Body(
        matter.Ground(),
        Transform(Vec3(0., 1.0, 0.)),
        aBody,
        Transform());
    cable.addViaPoint(cableViaPoint1Body, Vec3{0.});

    // Second via point.
    MobilizedBody::Translation cableViaPoint2Body(
        matter.Ground(),
        Transform(Vec3(0., 1.0, 2.0)),
        aBody,
        Transform());
    cable.addViaPoint(cableViaPoint2Body, Vec3{0.});

    // Third via point.
    MobilizedBody::Translation cableViaPoint3Body(
        matter.Ground(),
        Transform(Vec3(-4.0, 1.0, 2.0)),
        aBody,
        Transform());
    cable.addViaPoint(cableViaPoint3Body, Vec3{0.});

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
    cable.calcLength(s);
    if (viz) {
        viz->report(s);
    }

    SimTK_ASSERT_ALWAYS(
        cable.getNumObstacles() == 0,
        "Invalid number of obstacles");

    SimTK_ASSERT_ALWAYS(
        cable.getNumViaPoints() == 3,
        "Invalid number of via points");

    // With no obstacles, the smoothness should always be zero.
    SimTK_ASSERT1_ALWAYS(
        cable.getSmoothness(s) == 0.,
        "Test failed: Cable smoothness (=%e) must be zero",
        cable.getSmoothness(s));

    // Sum all straight line segment lengths.
    const Real lengthTolerance = 1e-5;
    const Real expectedTotalCableLength = 1. + 2. + 4. + std::sqrt(5.);
    const Real gotTotalCableLength = cable.calcLength(s);
    SimTK_ASSERT3_ALWAYS(
        std::abs(expectedTotalCableLength - gotTotalCableLength) <
            lengthTolerance,
        "Expected cable length (=%f) does not match computed cable length "
        "(=%f), error = %e",
        expectedTotalCableLength,
        gotTotalCableLength,
        expectedTotalCableLength - gotTotalCableLength);

    // We have no obstacles, so this should pass since all internal tests will
    // be skipped.
    CableSubsystemTestHelper().testCurrentPath(s, cables, std::cout);

    Vector_<SpatialVec> expectedBodyForcesInG(6, {{0., 0., 0.}, {0., 0., 0.}});
    // Ground body.
    expectedBodyForcesInG[0] = (SpatialVec{{0., 0., 0}, {0., 0., 0.}});
    // Cable origin body.
    expectedBodyForcesInG[1] = (SpatialVec{{0., 0., 0.}, {0., 1., 0.}});
    // Cable termination body.
    expectedBodyForcesInG[2] =
        (SpatialVec{{0., 0., 0.}, {0., 1./std::sqrt(5.), 2./std::sqrt(5.)}});
    // First via point.
    expectedBodyForcesInG[3] = (SpatialVec{{0., 0., 0.}, {0., -1., 1.}});
    // Second via point.
    expectedBodyForcesInG[4] = (SpatialVec{{0., 0., 0.}, {-1., 0, -1.}});
    // Third via point.
    expectedBodyForcesInG[5] =
        (SpatialVec{{0., 0., 0.}, {1., -1./std::sqrt(5.), -2./std::sqrt(5.)}});

    const Real tension = Pi;
    Vector_<SpatialVec> gotBodyForcesInG(6, {{0., 0., 0.}, {0., 0., 0.}});
    cable.applyBodyForces(s, tension, gotBodyForcesInG);

    for (int i = 0; i < gotBodyForcesInG.size(); ++i) {
        const Real tolerance = 1e-5;
        SimTK_ASSERT1_ALWAYS(
            (expectedBodyForcesInG[i] * tension - gotBodyForcesInG[i]).norm() <
                tolerance,
            "Unexpected force applied to body %i",
            i);
    }

    // Check tangent directions.
    // Origin.
    UnitVec3 originTangent = cable.calcOriginTangentDirection(s);
    SimTK_TEST_EQ(originTangent, UnitVec3(0., 1., 0.));

    // First via point.
    CableSpanViaPointIndex viaPoint1Index(0);
    UnitVec3 viaPoint1InTangent =
        cable.calcViaPointIncomingTangentDirection(s, viaPoint1Index);
    SimTK_TEST_EQ(viaPoint1InTangent, UnitVec3(0., 1., 0.));
    UnitVec3 viaPoint1OutTangent =
        cable.calcViaPointOutgoingTangentDirection(s, viaPoint1Index);
    SimTK_TEST_EQ(viaPoint1OutTangent, UnitVec3(0., 0., 1.));

    // Second via point.
    CableSpanViaPointIndex viaPoint2Index(1);
    UnitVec3 viaPoint2InTangent =
        cable.calcViaPointIncomingTangentDirection(s, viaPoint2Index);
    SimTK_TEST_EQ(viaPoint2InTangent, UnitVec3(0., 0., 1.));
    UnitVec3 viaPoint2OutTangent =
        cable.calcViaPointOutgoingTangentDirection(s, viaPoint2Index);
    SimTK_TEST_EQ(viaPoint2OutTangent, UnitVec3(-1., 0., 0.));

    // Third via point.
    CableSpanViaPointIndex viaPoint3Index(2);
    UnitVec3 viaPoint3InTangent =
        cable.calcViaPointIncomingTangentDirection(s, viaPoint3Index);
    SimTK_TEST_EQ(viaPoint3InTangent, UnitVec3(-1., 0., 0.));
    UnitVec3 viaPoint3OutTangent =
        cable.calcViaPointOutgoingTangentDirection(s, viaPoint3Index);
    SimTK_TEST_EQ(viaPoint3OutTangent,
                  UnitVec3(0., -1./std::sqrt(5.), -2./std::sqrt(5.)));

    // Termination.
    UnitVec3 terminationTangent = cable.calcTerminationTangentDirection(s);
    SimTK_TEST_EQ(terminationTangent,
                  UnitVec3(0., -1./std::sqrt(5.), -2./std::sqrt(5.)));
}

/** Test computed cable path over all supported surfaces and a via point,
testing geodesics, Jacobians, and kinematics.

The cable wraps over (in order):
1. Torus,
2. Ellipsoid,
3. Via point,
4. Sphere,
5. Cylinder,
6. Bicubic patch,
7. Torus.

The flag assertCableLengthDerivative is used to switch between verifying that
the computed cable length derivative matches the change in length during
simulation, OR verifying that all computed geodesics and Jacobians are
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

     // Mobilizer for the path via point.
    MobilizedBody::Translation cableViaPointBody(
        matter.Ground(),
        Transform(Vec3(0., 0.9, 0.5)),
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

    // Add the first via point.
    cable.addViaPoint(cableViaPointBody, Vec3{0.});

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
                4. * sin(angle * 0.7),
                10. * sin(angle * 1.3)));
        cableTerminationBody.setU(
            s,
            Vec3(
                0.1 * cos(angle),
                4. * 0.7 * cos(angle * 0.7),
                10. * 1.3 * cos(angle * 1.3)));

        // Move the via point.
        cableViaPointBody.setQ(
            s,
            Vec3(0., 0.5 * cos(angle), 0.));
        cableViaPointBody.setU(
            s,
            Vec3(0., 0.5 * -sin(angle), 0.));

        // Compute the CableSpan's path.
        system.realize(s, Stage::Report);
        const Real cableLength = cable.calcLength(s);
        if (viz) {
            viz->report(s);
        }

        // Check that the geodesics and path error vector & Jacobian are
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
            "Test failed: Cable length (=%f) smaller than distance between end "
            "points (=%f)",
            cableLength,
            distanceBetweenEndPoints);

        // Make sure that we acutally solved the path up to tolerance.
        SimTK_ASSERT2_ALWAYS(
            cable.getSmoothness(s) <= cable.getSmoothnessTolerance(),
            "Test failed: Cable smoothness (=%e) must be smaller than set "
            "tolerance (=%e)",
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

    // Helper for creating a cable with a single obstacle at a certain offset
    // location.
    auto createCable =
        [&](const std::function<ContactGeometry*()>& createSurface,
            const Transform& X_BS,
            const Transform& sceneOffset,
            bool useWarmStart            = true,
            const Vec3& contactPointHint = Vec3(NaN))
    {
        // Mobilizer for path origin.
        MobilizedBody::Translation cableOriginBody(
            matter.Ground(),
            sceneOffset.shiftFrameStationToBase(Vec3(-2., 0., 0.)),
            aBody,
            Transform());

        // Mobilizer for path termination.
        MobilizedBody::Translation cableTerminationBody(
            matter.Ground(),
            sceneOffset.shiftFrameStationToBase(Vec3(2., 0., 0.)),
            aBody,
            Transform());

        CableSpan cable(
            cables,
            cableOriginBody,
            Vec3{0.},
            cableTerminationBody,
            Vec3{0.});
        cable.setCurveSegmentAccuracy(1e-12);
        cable.setSmoothnessTolerance(1e-6);
        cable.setUseWarmStart(useWarmStart);

        if (contactPointHint.isFinite()) {
            cable.addObstacle(
                matter.Ground(),
                sceneOffset.compose(X_BS),
                std::shared_ptr<ContactGeometry>(createSurface()),
                contactPointHint);
        } else {
            cable.addObstacle(
                matter.Ground(),
                sceneOffset.compose(X_BS),
                std::shared_ptr<ContactGeometry>(createSurface()));
        }
    };

    // Create a cable with a torus obstacle.
    createCable(
        [&]() { return new ContactGeometry::Torus(2., 0.25); },
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{0., 1.75, 0.}),
        Transform(Vec3(0., 2., 0.)));

    // Create a cable with an ellipsoid obstacle.
    createCable(
        [&]() { return new ContactGeometry::Ellipsoid({1., 0.5, 0.75}); },
        Transform(Vec3{0., -0.5, 0.}),
        Transform(Vec3(0., 0., 0.)));

    // Create a cable with a sphere obstacle.
    createCable(
        [&]() { return new ContactGeometry::Sphere(1.5); },
        Transform(Vec3{0., -1.5, 0.}),
        Transform(Vec3(0., -2., 0.)));

    // Create a cable with a cylinder obstacle.
    createCable(
        [&]() { return new ContactGeometry::Cylinder(2.); },
        Transform(Vec3{0., -2., 0.}),
        Transform(Vec3(0., -6., 0.)));

    // Create a cable with a sphere obstacle, this time using a contact point
    // hint and no warm start. Without a warm start the path is recomputed from
    // the contact point hint at every realization, so liftoff is detected anew
    // each time the cable is out of contact. This is the configuration in which
    // a lifted curve segment is most likely to report a tracking point that is
    // still on the obstacle's surface rather than on the cable.
    createCable(
        [&]() { return new ContactGeometry::Sphere(1.5); },
        Transform(Vec3{0., -1.5, 0.}),
        Transform(Vec3(0., -10., 0.)),
        false,
        Vec3(0., 1.5, 0.));

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
        cables.getNumCables() == 5,
        "Unexpected number of cables (=%i)",
        cables.getNumCables());

    for (Real angle = 1e-2; angle < 4. * Pi; angle += 0.02) {
        // Move the cable end points.
        const Real yCoord = 0.1 * sin(angle);
        for (CableSpanIndex cableIx(0); cableIx < cables.getNumCables();
             ++cableIx) {
            matter
                .getMobilizedBody(cables.getCable(cableIx).getOriginBodyIndex())
                .setQToFitTranslation(s, Vec3(0., yCoord, 0.));
            matter
                .getMobilizedBody(
                    cables.getCable(cableIx).getTerminationBodyIndex())
                .setQToFitTranslation(s, Vec3(0., yCoord, 0.));
        }

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
                "Cable %i expected contact status (=%i) does not match "
                "computed contact status (=%i) at yCoord = %f",
                cableIx,
                expectedContactStatus,
                gotContactStatus,
                angle);

            // The decorative path points must be finite, whether or not the
            // obstacle is in contact with the cable.
            int numPoints = 0;
            auto assertPointIsFinite = [&](Vec3 point_G) {
                SimTK_ASSERT2_ALWAYS(point_G.isFinite(),
                    "Cable %i produced a non-finite decorative path point at "
                    "yCoord = %f", cableIx, yCoord);
                ++numPoints;
            };
            cables.getCable(cableIx).calcDecorativePathPoints(s,
                assertPointIsFinite);

            // The number of resampled points should remain the same, whether
            // or not the obstacle is in contact with the cable.
            constexpr int numSamples = 5;
            numPoints = 0;
            Array_<Vec3> resampledPoints;
            cables.getCable(cableIx).calcResampledDecorativePathPoints(s,
                numSamples, [&](Vec3 point_G) {
                    assertPointIsFinite(point_G);
                    resampledPoints.push_back(point_G);
                });
            SimTK_ASSERT3_ALWAYS(numPoints == numSamples + 2,
                "Cable %i expected %i resampled decorative path points, but "
                "got %i", cableIx, numSamples + 2, numPoints);

            // Each cable spans a single obstacle, so when that obstacle is out
            // of contact the cable is straight: every resampled point must lie
            // on the line connecting the origin and termination point.
            if (!gotContactStatus) {
                const Vec3 origin_G      = resampledPoints.front();
                const Vec3 line_G        = resampledPoints.back() - origin_G;
                const Real lineLengthSqr = line_G.normSqr();
                for (const Vec3& point_G : resampledPoints) {
                    const Real t =
                        dot(point_G - origin_G, line_G) / lineLengthSqr;
                    const Real distanceFromCable =
                        (point_G - (origin_G + t * line_G)).norm();
                    SimTK_ASSERT3_ALWAYS(distanceFromCable < 1e-9,
                        "Cable %i resampled decorative path point lies %g off "
                        "the cable at yCoord = %f",
                        cableIx, distanceFromCable, yCoord);
                }
            }
        }

        if (viz) {
            viz->report(s);
        }
    }
}

/** CableSpan path solver optimum test.

Span a cable over a single spherical obstacle, and check the solver result for
different algorithms:

- Algorithm::Scholz2015 will converge to the longest possible path.
- Algorithm::MinimumLength will converge to the shortest possible path. **/
void testSolverOptimum()
{
    const bool show = false;

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // Sphere obstacle radius.
    const Real radius = 1.;

    // The angle the cable makes w.r.t. the straight line connecting the
    // attachment points.
    const Real angle = 45. / 180. * Pi;

    // A dummy body.
    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    MobilizedBody::Free aMovingBody(
        matter.Ground(),
        Vec3(0.),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), radius, 0.},
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), -radius, 0.});

    // Add sphere obstacle, and choose the initial contact point near the
    // longest possible path.
    cable.addObstacle(
        aMovingBody,
        Transform(),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(radius)),
        Vec3(-radius, 0., 1e-1));

    // Configure the solver.
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-8);
    cable.setSolverMaxIterations(100);

    system.setUseUniformBackground(true); // no ground plane in display
    std::unique_ptr<Visualizer> viz(show ? new Visualizer(system) : nullptr);

    if (viz) {
        viz->setShowFrameNumber(true);
        viz->addDecorationGenerator(new CableDecorator(system, cable));
    }

    // Initialize the system and state.
    system.realizeTopology();

    // Scholz2015 algorithm converges to the longest possible path as the
    // optimal path.
    {
        cable.setAlgorithm(CableSpanAlgorithm::Scholz2015);

        State s = system.getDefaultState();
        system.realize(s, Stage::Report);

        if (viz) {
            viz->report(s);
        }

        SimTK_ASSERT2_ALWAYS(
            cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
            "Path smoothness (%e) does not meet tolerance (%e)",
            cable.getSmoothness(s),
            cable.getSmoothnessTolerance());

        const Vec3 expected_p_GP = Vec3(0., 1., 0.);
        const Vec3 got_p_GP      = cable
                                  .calcCurveSegmentInitialFrenetFrame(
                                      s,
                                      CableSpanObstacleIndex(0))
                                  .p();

        const Vec3 expected_p_GQ = Vec3(0., -1., 0.);
        const Vec3 got_p_GQ =
            cable.calcCurveSegmentFinalFrenetFrame(s, CableSpanObstacleIndex(0))
                .p();

        SimTK_ASSERT_ALWAYS(
            (expected_p_GP - got_p_GP).norm() < 1e-6,
            "Scholz2015 algorithm: Curve segment position at initial contact "
            "point incorrect");
        SimTK_ASSERT_ALWAYS(
            (expected_p_GQ - got_p_GQ).norm() < 1e-6,
            "Scholz2015 algorithm: Curve segment position at final contact "
            "point incorrect");
    }

    // MinimumLength algorithm converges to the shortest possible path as the
    // optimal path.
    {
        cable.setAlgorithm(CableSpanAlgorithm::MinimumLength);

        State s = system.getDefaultState();
        system.realize(s, Stage::Report);

        if (viz) {
            viz->report(s);
        }

        SimTK_ASSERT2_ALWAYS(
            cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
            "Path smoothness (%e) does not meet tolerance (%e)",
            cable.getSmoothness(s),
            cable.getSmoothnessTolerance());

        const Vec3 expected_p_GP =
            Vec3(cos(angle) * radius, sin(angle) * radius, 0.);
        const Vec3 got_p_GP = cable
                                  .calcCurveSegmentInitialFrenetFrame(
                                      s,
                                      CableSpanObstacleIndex(0))
                                  .p();

        const Vec3 expected_p_GQ =
            Vec3(cos(angle) * radius, -sin(angle) * radius, 0.);
        const Vec3 got_p_GQ =
            cable.calcCurveSegmentFinalFrenetFrame(s, CableSpanObstacleIndex(0))
                .p();

        SimTK_ASSERT_ALWAYS(
            (expected_p_GP - got_p_GP).norm() < 1e-6,
            "Scholz2015 algorithm: Curve segment position at initial contact "
            "point incorrect");
        SimTK_ASSERT_ALWAYS(
            (expected_p_GQ - got_p_GQ).norm() < 1e-6,
            "Scholz2015 algorithm: Curve segment position at final contact "
            "point incorrect");
    }
}

/** CableSpan robustness test during path initialization.

In this test the optimal path is far from the initial path. This tests the
robustness of the algorithm to the initial conditions. **/
void testRobustInitialPath()
{
    const bool show = false;

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // A dummy body.
    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    MobilizedBody::Free aMovingBody(
        matter.Ground(),
        Vec3(0.),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{-0.1, 0., 0.},
        matter.Ground(),
        Vec3{0.1, 0., 0.});

    // Add torus obstacle.
    cable.addObstacle(
        matter.Ground(),
        Transform(Rotation(0.5 * Pi, YAxis), Vec3{0., 4., 0.}),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Torus(2., 0.5)),
        {0., -0.1, 0.});

    // Configure the solver.
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-8);
    cable.setSolverMaxIterations(100);

    system.setUseUniformBackground(true); // no ground plane in display
    std::unique_ptr<Visualizer> viz(show ? new Visualizer(system) : nullptr);

    if (viz) {
        viz->setShowFrameNumber(true);
        viz->addDecorationGenerator(new CableDecorator(system, cable));
    }

    // Initialize the system and state.
    system.realizeTopology();

    // Try setting the algorithm to Algorithm::Scholz2015 to see how it fails
    // to initialize the path:
    // cable.setAlgorithm(CableSpanAlgorithm::Scholz2015);

    State s = system.getDefaultState();
    system.realize(s, Stage::Report);

    if (viz) {
        viz->report(s);
    }

    // Test that the solution was found.
    SimTK_ASSERT2_ALWAYS(
        cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
        "Path smoothness (%e) does not meet tolerance (%e)",
        cable.getSmoothness(s),
        cable.getSmoothnessTolerance());

    // Test the cable length.
    const Real expectedLength = 5.6513;
    SimTK_ASSERT2_ALWAYS(
        std::abs(cable.calcLength(s) - expectedLength) < 1e-3,
        "Cable length (%f) does match expected length (%f)",
        cable.calcLength(s),
        expectedLength);
}

/** CableSpan warm start test.

This tests that path solutions are dependent on previous path solutions when
warm starts are enabled, and only dependent on the contact hint when warm starts
are disabled. **/
void testUseWarmStart() {
    class PendulumSystem {
    public:
        PendulumSystem(bool useWarmStart, Real defaultAngle,
                Vec3 contactPointHint = Vec3(0.01, 0.2, 0.),
                int solverMaxIterations = 100,
                CableSpanAlgorithm algorithm =
                    CableSpanAlgorithm::MinimumLength) :
                m_matter(m_system), m_cables(m_system) {

            // Pendulum bodies.
            Body::Rigid pendulumBody(MassProperties(1., Vec3(0),
                    Inertia(1)));
            Body::Rigid headBody(MassProperties(1., Vec3(0), Inertia(1)));

            // Base body is welded to ground.
            MobilizedBody::Weld base(
                m_matter.Ground(),
                Transform(Vec3(0., 1.5, 0.)),
                pendulumBody,
                Transform());

            // Pendulum body is mobilized by a pin joint.
            m_pendulum = MobilizedBody::Pin(
                base,
                Transform(Vec3(0.)),
                headBody,
                Transform(Vec3(0., 1., 0.)));
            m_pendulum.setDefaultAngle(defaultAngle);

            // Add a cable from the head origin to the base origin, wrapping
            // over an ellipsoid obstacle.
            m_cable = CableSpan(
                m_cables,
                m_pendulum,
                Vec3(0.),
                base,
                Vec3(0.));
            // Add the ellipsoid obstacle with a contact hint biased in the
            // positive x-direction so that, when finding a path solution from
            // the contact hint, the cable will wrap over the side of the
            // ellipsoid in the positive x-direction.
            m_cable.addObstacle(
                base,
                Transform(Rotation(), Vec3(0, -0.5, 0)),
                std::shared_ptr<ContactGeometry>(
                    new ContactGeometry::Ellipsoid(Vec3(0.1, 0.1, 0.5))),
                contactPointHint);

            // Configure the solver.
            m_cable.setUseWarmStart(useWarmStart);
            m_cable.setCurveSegmentAccuracy(1e-12);
            m_cable.setSmoothnessTolerance(1e-8);
            m_cable.setSolverMaxIterations(solverMaxIterations);
            m_cable.setAlgorithm(algorithm);
        }

        State getDefaultState() {
            m_system.realizeTopology();
            return m_system.getDefaultState();
        }

        Real getCableLength(const State& s) {
            return m_cable.calcLength(s);
        }

        Vec3 getInitialCurvePoint(const State& s) {
            CableSpanObstacleIndex obsIx(0);
            return m_cable.calcCurveSegmentInitialFrenetFrame(s, obsIx).p();
        }

        void setPendulumAngle(State& s, Real angle) {
            m_pendulum.setAngle(s, angle);
        }

        void realizePosition(State& s) {
            m_system.realize(s, Stage::Position);
        }

        void assertConvergedForwardPath(const State& s) {
            SimTK_ASSERT_ALWAYS(
                m_cable.getSmoothness(s) < m_cable.getSmoothnessTolerance(),
                "The pendulum cable path should converge.");
            const Transform frame = m_cable.calcCurveSegmentFinalFrenetFrame(
                    s, CableSpanObstacleIndex(0));
            const Vec3 outgoingDirection = Vec3(0., 1.5, 0.) - frame.p();
            SimTK_ASSERT_ALWAYS(dot(frame.x(), outgoingDirection) > 0.,
                "The final wrap tangent should point toward the cable "
                "endpoint.");
            const Transform initialFrame =
                m_cable.calcCurveSegmentInitialFrenetFrame(
                    s, CableSpanObstacleIndex(0));
            const Vec3 incomingDirection =
                initialFrame.p() - m_pendulum.getBodyTransform(s).p();
            SimTK_ASSERT_ALWAYS(dot(initialFrame.x(), incomingDirection) > 0.,
                "The initial wrap tangent should point away from the cable "
                "origin.");
        }

        void setContactPointHint(const Vec3& hint) {
            m_cable.setObstacleContactPointHint(CableSpanObstacleIndex(0),
                hint);
        }

    private:
        MultibodySystem m_system;
        SimbodyMatterSubsystem m_matter;
        CableSubsystem m_cables;
        CableSpan m_cable;
        MobilizedBody::Pin m_pendulum;
    };


    // Given option to use warm start and a default pendulum angle, construct a
    // pendulum system with a wrap ellipsoid and calculate the initial curve
    // point. We first realize the system to the position stage using the
    // default angle, and then set the angle so that pendulum body to lies
    // directly benath the ellipsoid obstacle. The initial curve point is then
    // calculated. If warm start is enabled, the initial curve point should
    // depend on the choice of the default angle we dictates the initial path
    // solution. If warm start is disabled, the initial curve point should be
    // the same regardless of the default angle since the path solver should
    // always compute the path from the same initial point (i.e., the contact
    // hint).
    auto calcInitialCurvePoint = [](bool useWarmStart, Real defaultAngle,
            Real hintZ, CableSpanAlgorithm algorithm)
            -> Vec3 {
        PendulumSystem pendulumSystem(useWarmStart, defaultAngle,
                Vec3(0.01, 0.2, hintZ), 100, algorithm);
        State s = pendulumSystem.getDefaultState();
        pendulumSystem.realizePosition(s);
        pendulumSystem.assertConvergedForwardPath(s);
        pendulumSystem.setPendulumAngle(s, 2.0*SimTK::Pi);
        pendulumSystem.realizePosition(s);
        pendulumSystem.assertConvergedForwardPath(s);
        const Real radius = 0.1;
        const Real expectedLength = 2. * std::sqrt(0.25 - radius * radius)
                + 2. * radius * std::asin(radius / 0.5);
        SimTK_ASSERT_ALWAYS(
            std::abs(pendulumSystem.getCableLength(s) - expectedLength) < 1e-6,
            "The pendulum cable should follow the short, forward-facing wrap.");
        return pendulumSystem.getInitialCurvePoint(s);
    };


    // Both algorithms must find the same paths (see testSolverOptimum for a
    // case where they differ by design).
    const Real tol = 1e-5;
    for (CableSpanAlgorithm algorithm :
            {CableSpanAlgorithm::MinimumLength,
             CableSpanAlgorithm::Scholz2015}) {
        // Warm start enabled. The path should have converged to similar
        // solutions, but on opposite sides of the ellipsoid obstacle, so the
        // initial wrap points should have equal magnitude but opposite sign in
        // the x-direction and have equal y-components. The z-components should
        // be zero.
        {
            auto firstPoint =
                calcInitialCurvePoint(true, 0.05*SimTK::Pi, 0., algorithm);
            auto secondPoint =
                calcInitialCurvePoint(true, 1.95*SimTK::Pi, 0., algorithm);
            SimTK_ASSERT_ALWAYS(std::abs(firstPoint[0] + secondPoint[0]) < tol,
                "Warm start enabled: x-components of the initial wrap points "
                "should have equal magnitude but opposite sign, but they do "
                "not.");
            SimTK_ASSERT_ALWAYS(std::abs(firstPoint[1] - secondPoint[1]) < tol,
                 "Warm start enabled: y-components of the initial wrap points "
                 "should be equal, but they are not.");
            SimTK_ASSERT_ALWAYS(std::abs(firstPoint[2]) < tol &&
                                std::abs(secondPoint[2]) < tol,
                "Warm start enabled: z-components of the initial wrap points "
                "should be zero, but they are not.");
        }

        // Warm start disabled. The curve points should be equal since the
        // path solver should be computing the path from the same initial point
        // (i.e., the contact hint) in both case.
        {
            auto firstPoint =
                calcInitialCurvePoint(false, 0.05*SimTK::Pi, 0., algorithm);
            auto secondPoint =
                calcInitialCurvePoint(false, 1.95*SimTK::Pi, 0., algorithm);
            SimTK_ASSERT_ALWAYS((firstPoint - secondPoint).norm() < tol,
                "Warm start disabled: initial wrap points do not match.");
        }

        // Tiny out-of-plane perturbations must not send the warm-started path
        // toward the unstable semicircular wrap with a reversed terminal
        // tangent.
        for (bool useWarmStart : {false, true}) {
            const Vec3 point = calcInitialCurvePoint(useWarmStart,
                    1.95*SimTK::Pi, 1e-8, algorithm);
            SimTK_ASSERT_ALWAYS(std::abs(point[2]) < tol,
                "A perturbed contact hint should converge to the planar wrap.");
        }
    }

    // With no solver iterations, the initial point should be the nearest
    // surface projection of the hint, not the geodesic integrator's projection.
    // This must also hold for a hint that is changed after realizeTopology().
    const Vec3 hint(0.01, 0.2, 0.1);
    const Vec3 changedHint(-0.01, 0.2, -0.1);
    const ContactGeometry::Ellipsoid ellipsoid(Vec3(0.1, 0.1, 0.5));
    auto calcProjectedHint_G = [&](const Vec3& contactPointHint)
    {
        return ellipsoid.projectDownhillToNearestPoint(contactPointHint) +
               Vec3(0., 1., 0.);
    };
    for (bool useWarmStart : {false, true})
    for (bool changeHint : {false, true}) {
        PendulumSystem pendulumSystem(useWarmStart, 0.05*SimTK::Pi, hint, 0);
        State s = pendulumSystem.getDefaultState();
        if (changeHint) {
            pendulumSystem.setContactPointHint(changedHint);
        }
        pendulumSystem.realizePosition(s);
        const Vec3 expectedPoint_G =
            calcProjectedHint_G(changeHint ? changedHint : hint);
        SimTK_ASSERT_ALWAYS(
            (pendulumSystem.getInitialCurvePoint(s) - expectedPoint_G).norm()
                < tol,
            "The initial wrap point should use the nearest surface "
            "projection of the contact hint.");
    }
}

void testWarmStartNormalRestartDirection()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);
    Body::Rigid body(MassProperties(1., Vec3(0.), Inertia(1.)));
    MobilizedBody::Free termination(
        matter.Ground(), Transform(), body, Transform());
    termination.setDefaultTranslation(Vec3(1., -2., 0.));
    CableSpan cable(cables, matter.Ground(), Vec3(-2., 1., 0.),
            termination, Vec3(0.));
    cable.addObstacle(matter.Ground(), Transform(),
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Ellipsoid(Vec3(1., 1., 2.))),
        Vec3(0., 1., 0.));
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-8);
    cable.setSolverMaxIterations(100);
    cable.setAlgorithm(CableSpanAlgorithm::MinimumLength);

    system.realizeTopology();
    State s = system.getDefaultState();
    system.realize(s, Stage::Position);
    const CableSpanObstacleIndex obsIx(0);
    SimTK_ASSERT_ALWAYS(
        (cable.calcCurveSegmentInitialFrenetFrame(s, obsIx).p()
            - Vec3(0., 1., 0.)).norm() < 1e-6,
        "The initial cable should start its wrap at the top of the ellipsoid.");
    SimTK_ASSERT_ALWAYS(
        (cable.calcCurveSegmentFinalFrenetFrame(s, obsIx).p()
            - Vec3(1., 0., 0.)).norm() < 1e-6,
        "The initial cable should end its wrap at the side of the ellipsoid.");

    // Reverse the outgoing tangent while making next-minus-prev parallel
    // to the surface normal at the initial wrap point. No part of the folded
    // curve heads toward the new termination point, and no zero-length curve
    // at the initial wrap point can be aligned with both straight line
    // segments. The fold cannot be undone, so the curve must be left as is,
    // and the path must be reported as not smooth.
    cable.setSolverMaxIterations(0);
    const Transform initialFrame =
        cable.calcCurveSegmentInitialFrenetFrame(s, obsIx);
    termination.setQToFitTranslation(s,
        Vec3(-2., 1., 0.) + 2. * initialFrame.y());
    system.realize(s, Stage::Position);
    SimTK_ASSERT_ALWAYS(std::isfinite(cable.calcLength(s)),
        "A normal restart direction must not produce an invalid cable path.");
    SimTK_ASSERT_ALWAYS(cable.isInContactWithObstacle(s, obsIx),
        "The folded curve should be left in contact with the obstacle.");
    const Transform restartedFrame =
        cable.calcCurveSegmentInitialFrenetFrame(s, obsIx);
    SimTK_ASSERT_ALWAYS(dot(restartedFrame.x(), Vec3(1., 0., 0.)) > 0.99,
        "A normal restart direction should retain the finite initial tangent.");
    SimTK_ASSERT_ALWAYS(
        cable.getSmoothness(s) > cable.getSmoothnessTolerance(),
        "A folded cable path should not be reported as smooth.");
}

/** CableSpan fold recovery test.

A cable is hooked under a cylinder: both end points lie above the cylinder, and
the cable wraps more than 180 degrees around its underside. This is a valid
(locally shortest) path. Then the termination point is moved far in one step,
and the path is recomputed from the previous solution (warm start).

The previous geodesic now wraps too far around the cylinder, so that the
straight line leaving it runs backward along the curve. Such a fold has zero
normal and binormal path error, and the solver used to converge to it, reporting
it as a smooth path. The solver must instead find the hooked path through the
new termination point, for both algorithms.

All points lie in the plane z=0, so the expected path is a circular arc with
tangent lines that can be computed in closed form. **/
void testFoldRecoveryAfterLargeMotion()
{
    const Real radius = 1.;
    const Vec3 origin(-0.5, 2., 0.);
    const Vec3 initialTermination(0.5, 2., 0.);

    // Length of the path that leaves the origin along the left side of the
    // cylinder, wraps counterclockwise under it, and goes to the termination
    // point from its right side.
    auto calcHookedLength = [&](const Vec3& termination) -> Real
    {
        auto calcTangentPointAngle = [&](const Vec3& p, Real side) -> Real
        {
            return std::atan2(p[1], p[0]) + side * std::acos(radius / p.norm());
        };
        const Real angleIn = calcTangentPointAngle(origin, 1.);
        const Real angleOut = calcTangentPointAngle(termination, -1.);
        Real wrapAngle = std::fmod(angleOut - angleIn, 2. * Pi);
        if (wrapAngle < 0.) {
            wrapAngle += 2. * Pi;
        }
        return std::sqrt(origin.normSqr() - radius * radius) +
               radius * wrapAngle +
               std::sqrt(termination.normSqr() - radius * radius);
    };

    for (CableSpanAlgorithm algorithm :
            {CableSpanAlgorithm::MinimumLength,
             CableSpanAlgorithm::Scholz2015})
    for (const Vec3& newTermination :
            {Vec3(2., 0., 0.), Vec3(1.5, -1.5, 0.), Vec3(0., -2.5, 0.)}) {
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);
        Body::Rigid body(MassProperties(1., Vec3(0.), Inertia(1.)));
        MobilizedBody::Translation termination(
            matter.Ground(), Transform(), body, Transform());
        termination.setDefaultQ(initialTermination);

        CableSpan cable(cables, matter.Ground(), origin, termination, Vec3(0.));
        // The contact point hint below the cylinder selects the hooked path.
        cable.addObstacle(matter.Ground(), Transform(),
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Cylinder(radius)),
            Vec3(0., -1.5, 0.));
        cable.setCurveSegmentAccuracy(1e-12);
        cable.setSmoothnessTolerance(1e-8);
        cable.setSolverMaxIterations(100);
        // Scholz2015 cannot grow the initial zero-length curve into the
        // hooked path, so always use MinimumLength to find the initial path.
        cable.setAlgorithm(CableSpanAlgorithm::MinimumLength);

        system.realizeTopology();
        State s = system.getDefaultState();
        const CableSpanObstacleIndex obsIx(0);

        auto assertHookedPath = [&](const Vec3& terminationPoint)
        {
            system.realize(s, Stage::Position);
            SimTK_ASSERT_ALWAYS(
                cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
                "The hooked cable path should converge.");
            SimTK_ASSERT_ALWAYS(cable.isInContactWithObstacle(s, obsIx),
                "The hooked cable should remain in contact with the cylinder.");
            const Transform X_GP =
                cable.calcCurveSegmentInitialFrenetFrame(s, obsIx);
            const Transform X_GQ =
                cable.calcCurveSegmentFinalFrenetFrame(s, obsIx);
            SimTK_ASSERT_ALWAYS(dot(X_GP.x(), X_GP.p() - origin) > 0.,
                "The initial wrap tangent should point away from the origin.");
            SimTK_ASSERT_ALWAYS(dot(X_GQ.x(), terminationPoint - X_GQ.p()) > 0.,
                "The final wrap tangent should point toward the termination.");
            SimTK_ASSERT_ALWAYS(
                std::abs(cable.calcLength(s) -
                         calcHookedLength(terminationPoint)) < 1e-6,
                "The cable should follow the hooked path.");
        };

        // The initial hooked path.
        assertHookedPath(initialTermination);

        // Move the termination point in one step, and recompute the path from
        // the previous solution, using the algorithm under test.
        cable.setAlgorithm(algorithm);
        termination.setQ(s, newTermination);
        assertHookedPath(newTermination);
    }
    std::cout << "PASSED TEST: testFoldRecoveryAfterLargeMotion" << std::endl;
}

/** CableSpan fold recovery test for a cable that winds around an obstacle
multiple times.

A cable is wound more than once around a cylinder, by moving its termination
point around the cylinder in small steps. Then the termination point is moved
back far in one step, and the path is recomputed from the previous solution
(warm start). The previous geodesic now wraps too far around the cylinder, and
folds back at its final contact point. Undoing that fold must only unwind the
folded part, and keep the remaining full windings.

All points lie in the plane z=0, and the termination point stays at the same
distance from the cylinder axis, so the cable length changes by the angle that
the termination point moves (for a unit radius). **/
void testFoldRecoveryKeepsWindings()
{
    const Vec3 origin(1., -2., 0.);
    const Real terminationRadius = 2.;
    auto calcTermination = [&](Real angle) -> Vec3
    {
        return terminationRadius * Vec3(std::cos(angle), std::sin(angle), 0.);
    };

    for (CableSpanAlgorithm algorithm :
            {CableSpanAlgorithm::MinimumLength,
             CableSpanAlgorithm::Scholz2015}) {
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);
        Body::Rigid body(MassProperties(1., Vec3(0.), Inertia(1.)));
        MobilizedBody::Translation termination(
            matter.Ground(), Transform(), body, Transform());
        Real angle = 0.4 * Pi;
        termination.setDefaultQ(calcTermination(angle));

        CableSpan cable(cables, matter.Ground(), origin, termination, Vec3(0.));
        cable.addObstacle(matter.Ground(), Transform(),
            std::shared_ptr<ContactGeometry>(new ContactGeometry::Cylinder(1.)),
            Vec3(1., 0., 0.));
        cable.setCurveSegmentAccuracy(1e-12);
        cable.setSmoothnessTolerance(1e-8);
        cable.setSolverMaxIterations(100);
        cable.setAlgorithm(CableSpanAlgorithm::MinimumLength);

        system.realizeTopology();
        State s = system.getDefaultState();
        system.realize(s, Stage::Position);

        // Wind the cable around the cylinder, counterclockwise, until it
        // wraps about 3.2 pi.
        while (angle < 3.5 * Pi) {
            angle += 5. / 180. * Pi;
            termination.setQ(s, calcTermination(angle));
            system.realize(s, Stage::Position);
            SimTK_ASSERT_ALWAYS(
                cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
                "The cable should converge while winding around the cylinder.");
        }
        const Real woundLength = cable.calcLength(s);

        // Unwind by 0.9 pi in one step. More than two full windings remain.
        const Real unwindAngle = 0.9 * Pi;
        cable.setAlgorithm(algorithm);
        termination.setQ(s, calcTermination(angle - unwindAngle));
        system.realize(s, Stage::Position);
        SimTK_ASSERT_ALWAYS(
            cable.getSmoothness(s) < cable.getSmoothnessTolerance(),
            "The unwound cable path should converge.");
        SimTK_ASSERT_ALWAYS(
            std::abs(cable.calcLength(s) - (woundLength - unwindAngle)) < 1e-6,
            "Undoing the fold should keep the remaining windings.");
    }
    std::cout << "PASSED TEST: testFoldRecoveryKeepsWindings" << std::endl;
}

/** CableSpan exact fold test.

The cable starts with a zero-length curve on a cylinder whose straight line
segments run exactly in opposite directions along the curve tangent. The path
error of this fold is exactly zero, and the fold cannot be undone by restarting
the curve. The solver must not compute an invalid path correction from it (this
is checked by an assertion in Debug builds), and must either lift the cable off
the cylinder or report the path as not smooth. **/
void testExactFold()
{
    for (CableSpanAlgorithm algorithm :
            {CableSpanAlgorithm::MinimumLength,
             CableSpanAlgorithm::Scholz2015}) {
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);
        Body::Rigid body(MassProperties(1., Vec3(0.), Inertia(1.)));
        MobilizedBody::Translation termination(
            matter.Ground(), Transform(), body, Transform());
        termination.setDefaultQ(Vec3(-1., 1., 0.));

        CableSpan cable(
            cables, matter.Ground(), Vec3(-2., 1., 0.), termination, Vec3(0.));
        cable.addObstacle(matter.Ground(), Transform(),
            std::shared_ptr<ContactGeometry>(new ContactGeometry::Cylinder(1.)),
            Vec3(0., 1., 0.));
        cable.setCurveSegmentAccuracy(1e-12);
        cable.setSmoothnessTolerance(1e-8);
        cable.setSolverMaxIterations(20);
        cable.setAlgorithm(algorithm);

        system.realizeTopology();
        State s = system.getDefaultState();
        system.realize(s, Stage::Position);

        const Real length = cable.calcLength(s);
        SimTK_ASSERT_ALWAYS(std::isfinite(length),
            "An exact fold must not produce an invalid cable path.");
        if (cable.getSmoothness(s) < cable.getSmoothnessTolerance()) {
            SimTK_ASSERT_ALWAYS(
                !cable.isInContactWithObstacle(s, CableSpanObstacleIndex(0)) &&
                    std::abs(length - 1.) < 1e-10,
                "A smooth path must be the straight line between the end "
                "points.");
        }
    }
    std::cout << "PASSED TEST: testExactFold" << std::endl;
}

int main()
{
    testSimpleCable();
    testViaPoints();
    testTouchdownAndLiftoff();
    testAllSurfaceKinds(false); // Test all geodesics and Jacobians.
    testAllSurfaceKinds(true);  // Test length derivative.
    testSolverOptimum();
    testRobustInitialPath();
    testUseWarmStart();
    testWarmStartNormalRestartDirection();
    testFoldRecoveryAfterLargeMotion();
    testFoldRecoveryKeepsWindings();
    testExactFold();
}
