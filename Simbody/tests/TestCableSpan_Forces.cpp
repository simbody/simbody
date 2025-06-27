/*-----------------------------------------------------------------------------
                Simbody(tm) Test: Cable Forces
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

using namespace SimTK;

/**
This file contains simple tests for checking computed forces and moments.
**/

// Helper function for the force assertion.
void assertForces(
    MultibodySystem& system,
    CableSpan& cable,
    const Vector_<SpatialVec>& expectedBodyForcesInG)
{
    cable.setCurveSegmentAccuracy(1e-12);
    cable.setSmoothnessTolerance(1e-8);

    // Initialize the system and state.
    system.realizeTopology();
    State s = system.getDefaultState();
    system.realize(s, Stage::Report);

    const int numBodies = system.getMatterSubsystem().getNumBodies();
    Vector_<SpatialVec> forces(numBodies, SpatialVec{Vec3{0.}, Vec3{0.}});
    cable.applyBodyForces(s, 1., forces);

    SimTK_ASSERT2_ALWAYS(
        expectedBodyForcesInG.size() == numBodies,
        "Size of expected forces (=%i) does not match number of bodies (=%i)",
        expectedBodyForcesInG.size(),
        numBodies);

    for (int i = 0; i < numBodies; ++i) {
        SimTK_ASSERT1_ALWAYS(
            (forces[i] - expectedBodyForcesInG[i]).norm() < 1e-6,
            "Test failed: force error = %e",
            (forces[i] - expectedBodyForcesInG[i]).norm());
    }
}

// A cable with both end points on the same body, and no obstacles, should
// result in zero net force.
void testCableForceOnSameBody()
{
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // Construct a new cable.
    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{0., 1., 0.},
        matter.Ground(),
        Vec3{1., 1., 0.});

    // The expected result.
    const Vector_<SpatialVec> expectedBodyForcesInG(
        1,
        SpatialVec(Vec3{0.}, Vec3{0.}));
    assertForces(system, cable, expectedBodyForcesInG);
}

// A cable spanned between two bodies, and no obstacles, should result in a
// force pulling them towards each other.
void testCableForceBetweenTwoBodies()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    MobilizedBody::Free cableOriginBody(
        matter.Ground(),
        Vec3(1., 0., 0.),
        aBody,
        Transform());

    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{0.},
        cableOriginBody,
        Vec3{0., 0., 0.});

    Vector_<SpatialVec> expectedBodyForcesInG(matter.getNumBodies());
    expectedBodyForcesInG[0] = SpatialVec(Vec3{0.}, Vec3{1., 0., 0.});
    expectedBodyForcesInG[1] = SpatialVec(Vec3{0.}, Vec3{-1., 0., 0.});
    assertForces(system, cable, expectedBodyForcesInG);
}

// Consider a cable spanned between two bodies and no obstacles: If the
// attachment point is not at the body origin, it should generate a moment on
// the body.
void testCableForceAndMomentBetweenTwoBodies()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    MobilizedBody::Free cableOriginBody(
        matter.Ground(),
        Vec3(1., -2., 0.),
        aBody,
        Transform());

    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{0.},
        cableOriginBody,
        Vec3{0., 2., 0.});

    Vector_<SpatialVec> expectedBodyForcesInG(matter.getNumBodies());
    expectedBodyForcesInG[0] = SpatialVec(Vec3{0.}, Vec3{1., 0., 0.});
    expectedBodyForcesInG[1] = SpatialVec(Vec3{0., 0., 2.}, Vec3{-1., 0., 0.});
    assertForces(system, cable, expectedBodyForcesInG);
}

// Imagine a sling-shot: A cable with both attachment points fixed in Ground,
// and a sphere attached to an obstacle pulling on the cable in the middle. The
// sphere is positioned in such a way that the cable makes a 45 degree angle
// w.r.t. the straight line connecting the attachment points.
void testCableForceOnSingleObstacle()
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
    MobilizedBody::Free obstacleBody(
        matter.Ground(),
        Vec3(0.),
        aBody,
        Transform());

    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), radius, 0.},
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), -radius, 0.});

    // Add sphere obstacle.
    cable.addObstacle(
        obstacleBody,
        Transform(),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(radius)),
        Vec3(radius, 0., 0.));

    Vector_<SpatialVec> forcesExpected(matter.getNumBodies());
    forcesExpected[0] = SpatialVec(Vec3{0.}, Vec3{2. * cos(angle), 0., 0.});
    forcesExpected[1] =
        SpatialVec(Vec3{0., 0., 0.}, Vec3{-2. * cos(angle), 0., 0.});

    assertForces(system, cable, forcesExpected);
}

// This case is the same as testCableForceOnSingleObstacle, except that there is 
// an offset between the obstacle surface frame and obstacle body, to verify the 
// generated torque.
void testCableForceAndMomentOnSingleObstacle()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));

    const Real radius  = 1.;
    const Real offsetY = 1.;
    const Real angle   = 45. / 180. * Pi;

    MobilizedBody::Free obstacleBody(
        matter.Ground(),
        Vec3(0., -offsetY, 0.),
        aBody,
        Transform());

    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), radius, 0.},
        matter.Ground(),
        Vec3{radius * tan(angle / 2.), -radius, 0.});

    cable.addObstacle(
        obstacleBody,
        Transform(Vec3(0., offsetY, 0.)),
        std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(radius)),
        Vec3(radius, 0., 0.));

    Vector_<SpatialVec> forcesExpected(matter.getNumBodies());
    forcesExpected[0] = SpatialVec(Vec3{0.}, Vec3{2. * cos(angle), 0., 0.});
    forcesExpected[1] = SpatialVec(
        Vec3{0., 0., 2. * sin(angle)},
        Vec3{-2. * cos(angle), 0., 0.});

    assertForces(system, cable, forcesExpected);
}

// Similar to testCableForceOnSingleObstacle, but with a via point pulled
// along the X-axis instead of an obstacle.
void testCableForceOnSingleViaPoint()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // The angle the cable makes w.r.t. the straight line connecting the
    // attachment points.
    const Real angle = 45. / 180. * Pi;

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    
    MobilizedBody::Free viaPointBody(
        matter.Ground(),
        Vec3(0.),
        aBody,
        Transform());

    // Construct a new cable.
    CableSpan cable(
        cables,
        matter.Ground(),
        Vec3{0., 1., 0.},
        matter.Ground(),
        Vec3{0., -1., 0.});

    // Add via point.
    cable.addViaPoint(viaPointBody, Vec3(tan(angle), 0., 0.));

    Vector_<SpatialVec> forcesExpected(matter.getNumBodies());
    forcesExpected[0] = SpatialVec(Vec3{0.}, Vec3{2. * sin(angle), 0., 0.});
    forcesExpected[1] = SpatialVec(Vec3{0.}, Vec3{-2. * sin(angle), 0., 0.});

    assertForces(system, cable, forcesExpected);
}

// Similar to testCableForceAndMomentOnSingleObstacle, but with a via point
// pulled along the X-axis instead of an obstacle.
void testCableForceAndMomentOnSingleViaPoint()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableSubsystem cables(system);

    // The angle the cable makes w.r.t. the straight line connecting the
    // attachment points.
    const Real angle = 45. / 180. * Pi;

    Body::Rigid aBody(MassProperties(1., Vec3(0), Inertia(1)));
    const Real offsetY = 1.;
    MobilizedBody::Free viaPointBody(
        matter.Ground(),
        Vec3(0., -offsetY, 0.),
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
    cable.addViaPoint(viaPointBody, Vec3(tan(angle), offsetY, 0.));

    Vector_<SpatialVec> forcesExpected(matter.getNumBodies());
    forcesExpected[0] = SpatialVec(Vec3{0.}, Vec3{2. * sin(angle), 0., 0.});
    forcesExpected[1] = SpatialVec(
        Vec3{0., 0., 2. * sin(angle)},
        Vec3{-2. * sin(angle), 0., 0.});

    assertForces(system, cable, forcesExpected);
}

int main()
{
    testCableForceOnSameBody();
    testCableForceBetweenTwoBodies();
    testCableForceAndMomentBetweenTwoBodies();
    testCableForceOnSingleObstacle();
    testCableForceAndMomentOnSingleObstacle();
    testCableForceOnSingleViaPoint();
    testCableForceAndMomentOnSingleViaPoint();
}
