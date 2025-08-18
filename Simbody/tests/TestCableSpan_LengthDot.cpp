/*-----------------------------------------------------------------------------
                Simbody(tm) Test: Cable Lengthening Speed
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
This file contains simple tests for checking computed path lengthening speeds.
**/

// Helper function for the force assertion.
void assertLengthDot(
    const std::string& testCase,
    MultibodySystem& system,
    CableSpan& cable,
    const std::function<void(State& s)>& setU,
    Real expectedLengthDot)
{
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-8);

    // Initialize the system and state.
    system.realizeTopology();
    State s = system.getDefaultState();
    setU(s);
    system.realize(s, Stage::Report);

    const Real gotLengthDot = cable.calcLengthDot(s);

    for (CableSpanObstacleIndex ix(0); ix < cable.getNumObstacles(); ++ix) {
        SimTK_ASSERT2_ALWAYS(
            cable.isInContactWithObstacle(s, ix),
            "%s failed: Cable not in contact with obstacle %i",
            testCase.c_str(),
            ix);
    }

    SimTK_ASSERT4_ALWAYS(
        std::abs(gotLengthDot - expectedLengthDot) < 1e-13,
        "%s failed: expected lengthDot (=%f) does not match computed "
        "lengthDot (=%f), with error %e",
        testCase.c_str(),
        expectedLengthDot,
        gotLengthDot,
        expectedLengthDot - gotLengthDot);
}

// A cable with both end points on the same body, and no obstacles, should
// result in zero lengthening speed.
void testCableLengthDotOnSameBody()
{
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    MobilizedBody::Free aMovingBody(
        matter.Ground(),
        Vec3(0., 0., 0.),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        aMovingBody,
        Vec3(-1., 0., 0.),
        aMovingBody,
        Vec3(0., 1., 0.));

    // Moving body velocity:
    const Vec6 u{
        // angular velocity
        1., 2., 3.,
        // linear velocity
        4., 5., 6.,
    };
    const Real expectedLengthDot = 0.;

    assertLengthDot(
        "testCableLengthDotOnSameBody",
        system,
        cable,
        [&](State& s) { aMovingBody.setU(s, u); },
        expectedLengthDot);
}

// A cable with both end points on the same body, and two obstacles on that
// body, should result in zero lengthening speed.
void testCableLengthDotOnSameBodyWitObstacles()
{
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    MobilizedBody::Free aMovingBody(
        matter.Ground(),
        Transform(),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        aMovingBody,
        Vec3(2., 0., 0.),
        aMovingBody,
        Vec3(-2., 0., 0.));

    cable.addObstacle(
        aMovingBody,
        Transform(Vec3(1., 1., 0.)),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(1.)),
        Vec3(1., 1., 0.));

    cable.addObstacle(
        aMovingBody,
        Transform(Vec3(-1., 1., 0.)),
        std::shared_ptr<ContactGeometry>(
            new ContactGeometry::Ellipsoid(Vec3{0.5, 0.75, 1.})),
        Vec3(-1., 1., 0.));

    const Vec6 u{
        1., 2., 3.,
        4., 5., 6.,
    };
    const Real expectedLengthDot = 0.;

    assertLengthDot(
        "testCableLengthDotOnSameBodyWitObstacles",
        system,
        cable,
        [&](State& s) { aMovingBody.setU(s, u); },
        expectedLengthDot);
}

// Consider a cable spanned between two bodies and no obstacles. Verify
// lengthening speed for different body velocities.
void testCableLengthDotBetweenTwoBodies()
{
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    MobilizedBody::Free aMovingBody(
        matter.Ground(),
        Vec3(1., -2., 0.),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3(0., 0., 0.),
        aMovingBody,
        Vec3(0., 2., 0.));

    // Translation only.
    {
        const Vec6 u{
            0., 0., 0.,
            4., 5., 6.,
        };
        // Since the cable direction is along the XAxis:
        const Real expectedLengthDot = u[3];

        assertLengthDot(
            "testCableLengthDotBetweenTwoBodies(Translation)",
            system,
            cable,
            [&](State& s) { aMovingBody.setU(s, u); },
            expectedLengthDot);
    }

    // Rotation only.
    {
        const Vec6 u{
            0., 1., 2.,
            0., 0., 0.,
        };
        // Cross the angular rate with the attachment point offset.
        const Real expectedLengthDot = -2. * 2.;

        assertLengthDot(
            "testCableLengthDotBetweenTwoBodies(Rotation)",
            system,
            cable,
            [&](State& s) { aMovingBody.setU(s, u); },
            expectedLengthDot);
    }

    // Both.
    {
        const Vec6 u{
            0., 1., 2.,
            4., 5., 6.,
        };
        const Real expectedLengthDot = 0.;

        assertLengthDot(
            "testCableLengthDotBetweenTwoBodies(Both)",
            system,
            cable,
            [&](State& s) { aMovingBody.setU(s, u); },
            expectedLengthDot);
    }
}

// Sling shot case: A cable with both attachment points fixed in Ground, and a
// sphere obstacle which is pulled along the X-axis. The obstacle causes the
// cable to make a 45 degree angle w.r.t. the straight line connecting the
// attachment points.
void testCableLengthDotWithSingleObstacle()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // Sphere obstacle radius.
    const Real radius = 1.;

    // The angle the cable makes w.r.t. the straight line connecting the
    // attachment points.
    const Real angle = 45. / 180. * Pi;

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

    // Add sphere obstacle.
    cable.addObstacle(
        aMovingBody,
        Transform(),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(radius)),
        Vec3(radius, 0., 0.));

    const Vec6 u{
        1., 2., 3.,
        4., 5., 6.,
    };
    // Due to symmetry only the x-component of the obstacle velocity
    // contributes to the lengthening speed, with the other components
    // canceling out.
    const Real expectedLengthDot = 4. * 2. * cos(angle);

    assertLengthDot(
        "testCableLengthDotWithSingleObstacle",
        system,
        cable,
        [&](State& s) { aMovingBody.setU(s, u); },
        expectedLengthDot);
}

// Similar to testCableLengthDotWithSingleObstacle, but with a via point pulled
// along the X-axis instead of an obstacle.
void testCableLengthDotWithSingleViaPoint()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // The angle the cable makes w.r.t. the straight line connecting the
    // attachment points.
    const Real angle = 45. / 180. * Pi;

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
        Vec3{0, 1., 0.},
        matter.Ground(),
        Vec3{0, -1., 0.});

    // Add via point.
    cable.addViaPoint(aMovingBody, Vec3(tan(angle), 0., 0.));

    const Vec6 u{
        1., 2., 3.,
        4., 5., 6.,
    };
    // Due to symmetry only the x-component of the via point velocity
    // contributes to the lengthening speed, with the other components
    // canceling out.
    const Real expectedLengthDot = 4. * 2. * sin(angle);

    assertLengthDot(
        "testCableLengthDotWithSingleViaPoint",
        system,
        cable,
        [&](State& s) { aMovingBody.setU(s, u); },
        expectedLengthDot);
}

int main()
{
    testCableLengthDotOnSameBody();
    testCableLengthDotOnSameBodyWitObstacles();
    testCableLengthDotBetweenTwoBodies();
    testCableLengthDotWithSingleObstacle();
    testCableLengthDotWithSingleViaPoint();
}
