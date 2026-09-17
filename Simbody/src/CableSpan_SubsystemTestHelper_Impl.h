#ifndef SimTK_SIMBODY_CABLE_SUBSYSTEM_TEST_HELPER_H_
#define SimTK_SIMBODY_CABLE_SUBSYSTEM_TEST_HELPER_H_

/*-----------------------------------------------------------------------------
                               Simbody(tm)
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

#include "simbody/internal/CableSpan.h"

namespace SimTK
{

/** This class is the implementation side of CableSubsystemTestHelper. **/
class CableSubsystemTestHelper::Impl final {
public:
    //--------------------------------------------------------------------------
    //  Interface Translation Function
    //--------------------------------------------------------------------------

    /** See CableSpan::testCurrentPath **/
    void testCurrentPath(
        const State& state,
        const CableSubsystem& subsystem,
        std::ostream& testReport) const;

    //--------------------------------------------------------------------------
    // Configuration Parameters
    //--------------------------------------------------------------------------
    // These data structures contain the parameters for configuring the tests.

    /** Parameters for configuring the geodesic test. **/
    struct GeodesicTestParameters {
        /** Integrator accuracy of the unconstrained geodesic. **/
        Real accuracy = 1e-13;
        /** Tolerance used to assert constraint violation after projecting the
        curve segment's initial conditions to the surface. **/
        Real constraintTolerance = 1e-9;
        /** Tolerance used to assert constraint violation of the unconstrained
        geodesic state at the final arc length. **/
        Real driftTolerance = 1e-8;
        /** Tolerance used to assert equality between the CableSpan's curve
        segment and the unconstrained geodesic. **/
        Real assertionTolerance = 1e-7;
    };

    /** Parameters for configuring the jacobian perturbation test. **/
    struct PerturbationTestParameters {
        /** Perturbation value used. **/
        Real perturbation = 1e-5;
        /** Tolerance used to assert that the jacobian correctly predicts the
        perturbation effect. This is a percentage: e.g. tolerance = 1 would
        mean that the jacobian has a prediction error less or equal to 1
        percent for the given perturbation. **/
        Real tolerance = 0.1;
    };

private:
    //--------------------------------------------------------------------------
    // Tests
    //--------------------------------------------------------------------------

    /** This function asserts that the internally computed path error jacobian
    is correct using a perturbation test.

    The perturbation test is done as follows:
    1. Compute the path error vector and jacobian before perturbation.
    2. Take the perturbation as a small nonzero NaturalGeodesicCorrection, for a
        single obstacle, and predict the effect on the path error vector using
        the jacobian.
    3. Apply the NaturalGeodesicCorrection (i.e. the perturbation) to the path,
        and evaluate the path error vector again.
    4. The test passes if the predicted path error vector obtained from the
        jacobian, matches the evaluated path error vector after the
        perturbation, up to tolerance.

    The above test is repeated for each of the four NaturalGeodesicCorrection
    elements, and for each obstacle that comes into contact with the cable.

    Note that the jacobian is a local predictor of the effect of the
    NaturalGeodesicCorrection on the path error vector. This means that
    switching behavior due to touchdown and liftoff would result in an incorrect
    prediction. Whenever a switching behavior is detected, the test is skipped.

    Interesting test results are written to the testReport.
    An exception is thrown when the test fails. **/
    void runPerturbationTest(
        const State& state,
        const CableSubsystem& subsystem,
        std::ostream& testReport) const;

    /** This function asserts that each curve segment of each CableSpan in the
    CableSubsystem is a correct geodesic.

    This test will verify correctness of the ContactGeometry's computed:
    - normal curvature
    - geodesic torsion
    - surface constraint gradient direction
    and verify that any curve segment of the CableSpan's path is a geodeisc.
    It is assumed that ContactGeometry::calcSurfaceValue is correct.

    Verifying that the segment is a geodesic is done by computing an
    unconstrained geodesic starting from the same initial conditions. The
    unconstrained geodesic uses the implicit geodesic equations to compute the
    frenet frame at the final arc length, but does not enforce any constraints
    during the numerical integration.

    The unconstrained geodesic uses the following state:
    - x: The point, or position in ground
    - t: The tangent vector in ground
    - n: The normal vector in ground
    - b: The binormal vector in ground

    Consider the differential equations for integrating the position of the
    point on the geodesic:

    dx = t
    dt = f(x) = kn(x) * g(x) / |g(x)|

    where kn(x) is the normal curvature, and g(x) is the gradient of the surface
    constraint. Clearly the dynamics of the point and tangent form an isolated
    system of equations along the geodesic. If we integrate this system to the
    final arc length, without enforcing any constraints (i.e. no surface
    projection, tangent norm, tangent plane, etc), then any bugs in either kn(x)
    or g(x) will cause the constraints to be violated. These constraint
    violations are easily evaluated:
    1. The point should not drift away from the surface.
    2. The tangent vector should have unit norm.
    3. The tangent vector should be perpendicular to the surface normal.
    So, we can test the normal curvature and gradient direction using the
    position and tangent of the final frenet frame, of the unconstrained
    geodesic.

    Next, consider the differential equation for the binormal vector:

    db = tau_g(x, t) * g(x) / |g(x)|

    where tau_g(x, t) is the geodesic torsion. Again we can evaluate the
    binormal after integrating to the final arc length: It should have unit
    length, and be perpendicular to the tangent and surface normal. Any bugs in
    computing the geodesic torsion will result in a strange direction of the
    binormal at the final arc length.

    Finally, we have the differential equation for the normal vector:

    dn = -tau_g * b - kn * t

    Integrating this differential equation should result in a vector that
    matches computing the surface normal at the final arc length.

    Given the initial conditions, the arc length, the normal curvature and
    geodesic torsion: the geodesic is completely determined. We assert that the
    unconstrained geodesic gives the same final frenet frame as the curve
    segment over the obstacle, and we can be sure that this curve segment is
    indeed a geodesic.

    Interesting test results are written to the testReport.
    An exception is thrown when the test fails. **/
    void runGeodesicTest(
        const State& state,
        const CableSubsystem& subsystem,
        std::ostream& testReport) const;

    //--------------------------------------------------------------------------
    //  Data: Configuration parameters
    //--------------------------------------------------------------------------
    // NOTE Currently these are not exposed, so they are effectively fixed to
    // their default values. But should this ever change, it helps that the
    // current structure is in place.

    /* Configuration parameters for the jacobian perturbation test. */
    PerturbationTestParameters m_perturbationTestParameters;

    /* Configuration parameters for the geodesic test. */
    GeodesicTestParameters m_integratorTestParameters;
};

} // namespace SimTK

#endif
