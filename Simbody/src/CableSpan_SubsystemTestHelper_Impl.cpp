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

/* This file contains the implementation of the geodesic test for the curve
segments of a CableSpan. The implementation of the jacobian perturbation test
can be found in the CableSpan::Impl source file, because the jacobian
calculations are not exposed through the public interface.

We will use the following notation:
- P subscript: First contact point on obstacle
- Q subscript: Final contact point on obstacle
- G subscript: Ground frame
- S subscript: Surface frame
- X_GP, X_GQ: Frenet frame at first and final contact point w.r.t. ground.
- X_GS: Position and orientation of surface frame w.r.t. ground.
- x: Position at point on geodesic.
- t: Tangent at point on geodesic.
- n: Normal at point on geodesic.
- b: Binormal at point on geodesic. */

#include "CableSpan_SubsystemTestHelper_Impl.h"

#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simmath/internal/ContactGeometry.h"

using namespace SimTK;

namespace
{
using IntegratorTestParameters =
    CableSubsystemTestHelper::Impl::GeodesicTestParameters;
using PerturbationTestParameters =
    CableSubsystemTestHelper::Impl::PerturbationTestParameters;
} // namespace

//==============================================================================
//              Unconstrained Geodesic for the Geodesic Test
//==============================================================================
/* The following section defines the unconstrained geodesic which is used to
test the normal curvature, geodesic torsion, and surface gradient. Furthermore,
we can use the unconstrained geodesic to verify that the curve segments of a
CableSpan are indeed geodesics. */
namespace
{

/* This struct contains the state of the unconstrained geodesic. */
struct UnconstrainedGeodesicState {
    // Position of the current point on the geodeisc.
    Vec3 x{NaN};
    // The (unconstrained) Frenet frame orientation at the current point on the
    // geodesic, with the columns containing the tangent, normal and binormal
    // directions (in that order). Note that this orientation should be a
    // Rotation matrix, but since it is unconstrained, it is just a Mat33.
    Mat33 R{NaN}; // = [t_G, n_G, b_G]
};

/* This helper converts an UnconstrainedGeodesicState to a Frenet frame. */
Transform calcFrenetFrame(const UnconstrainedGeodesicState& q)
{
    const UnitVec3 n(q.R.col(1));
    const Vec3 t = q.R.col(0);
    return {Rotation().setRotationFromTwoAxes(n, YAxis, t, XAxis), q.x};
}

/* Helper function for projecting the approximate initial Frenet frame
X_GP_approx to the surface. It is important when integrating the
UnconstrainedGeodesicState that we do start on the surface. */
Transform projectFirstContactFrameToSurface(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    const Transform& X_GP_approx)
{
    const Vec3 x_S = geometry.projectDownhillToNearestPoint(
        X_GS.shiftBaseStationToFrame(X_GP_approx.p()));
    const UnitVec3 n_S = geometry.calcSurfaceUnitNormal(x_S);

    const Vec3 x_G(X_GS.shiftFrameStationToBase(x_S));
    const UnitVec3 n_G(X_GS.xformFrameVecToBase(n_S));

    return {
        Rotation().setRotationFromTwoAxes(
            n_G,
            YAxis,
            X_GP_approx.R().getAxisUnitVec(XAxis),
            XAxis),
        x_G};
}

/* Helper function for computing a CableSpan's obstacle's surface to ground
transform. */
Transform getXformSurfaceToGround(
    const State& state,
    const CableSubsystem& subsystem,
    CableSpanIndex cableIx,
    CableSpanObstacleIndex obsIx)
{
    const CableSpan& cable = subsystem.getCable(cableIx);
    const Mobod& body =
        subsystem.getMultibodySystem().getMatterSubsystem().getMobilizedBody(
            cable.getObstacleMobilizedBodyIndex(obsIx));
    return body.getBodyTransform(state).compose(
        cable.getObstacleTransformSurfaceToBody(obsIx));
}

/* This is the unconstrained geodesic (system) which has the
UnconstrainedGeodesicState as state, and uses the RungeKuttaMersonIntegrator to
integrate the state to the final arc length.

The elements in the State vector are ordered as:

S = [x_G, t_G, n_G, b_G]

The goal is to verify correctness of the CableSpan's curve segment over an
obstacle: Therefore the constructor takes the obstacle geometry, the surface
offset frame, and the curve segment's initial Frenet frame. With this
information we should be able to compute the geodesic that matches the curve
segment.
*/
class UnconstrainedGeodesic : public System {
public:
    /* Constructor that takes all information to compute the geodesic that
    matches the CableSpan's curve segment. */
    UnconstrainedGeodesic(
        const ContactGeometry& geometry,
        const Transform& X_GS,
        const Transform& X_GP);

    /* Construct the UnconstrainedGeodesic that will match the curve segment of
    the given cable over the given obstacle */
    UnconstrainedGeodesic(
        const State& state,
        const CableSubsystem& subsystem,
        CableSpanIndex cableIx,
        CableSpanObstacleIndex obsIx) :
        UnconstrainedGeodesic(
            subsystem.getCable(cableIx).getObstacleContactGeometry(obsIx),
            getXformSurfaceToGround(state, subsystem, cableIx, obsIx),
            subsystem.getCable(cableIx).calcCurveSegmentInitialFrenetFrame(
                state,
                obsIx))
    {}

    /* Helper function for getting the UnconstrainedGeodesicState from the
    State vector. */
    static UnconstrainedGeodesicState getGeodesicState(const State& state)
    {
        UnconstrainedGeodesicState y;
        int ix = -1;

        for (int row = 0; row < 3; ++row) {
            y.x[row] = state.getZ()[++ix];
        }

        for (int col = 0; col < 3; ++col) {
            for (int row = 0; row < 3; ++row) {
                y.R(row, col) = state.getZ()[++ix];
            }
        }

        return y;
    }

    class Guts : public System::Guts {
        friend UnconstrainedGeodesic;

        explicit Guts(
            const ContactGeometry& geometry,
            const Transform& X_GS,
            const Transform& X_GP) :
            m_geometry(geometry), m_X_GS(X_GS),
            m_X_GP(projectFirstContactFrameToSurface(geometry, m_X_GS, X_GP))
        {}

        Guts* cloneImpl() const override
        {
            return new Guts(*this);
        }

        /* During realizeTopology() we allocate the State for holding the
        UnconstrainedGeodesicState. */
        int realizeTopologyImpl(State& state) const override
        {
            Vector zInit(3 + 9, 0.);

            int ix = -1;

            for (int row = 0; row < 3; ++row) {
                zInit[++ix] = m_X_GP.p()[row];
            }

            for (int col = 0; col < 3; ++col) {
                for (int row = 0; row < 3; ++row) {
                    zInit[++ix] = m_X_GP.R().asMat33()(row, col);
                }
            }

            state.allocateZ(SubsystemIndex(0), zInit);
            return 0;
        }

        // Calculates the State (=UnconstrainedGeodesicState) derivative.
        int realizeAccelerationImpl(const State& state) const override
        {
            UnconstrainedGeodesicState y = getGeodesicState(state);

            const Vec3 t_G = y.R.col(0);

            const Vec3 x_S = m_X_GS.shiftBaseStationToFrame(y.x);
            const UnitVec3 t_S(m_X_GS.xformBaseVecToFrame(t_G));

            // NOTE do not use the normal from the state, but compute it from
            // the point, as explained in
            // CableSubsystemTestHelper::Impl::runGeodesicTest.
            const UnitVec3 n_S = m_geometry.calcSurfaceUnitNormal(x_S);

            const UnitVec3 n_G(m_X_GS.xformFrameVecToBase(n_S));
            const UnitVec3 b_G(t_G % n_G);

            // NOTE Simbody uses curvature with sign flipped from convention.
            const Real kn =
                -m_geometry.calcSurfaceCurvatureInDirection(x_S, t_S);
            const Real tau_g =
                -m_geometry.calcSurfaceTorsionInDirection(x_S, t_S);

            int ix = -1;

            // xDot
            state.updZDot()[++ix] = t_G[0];
            state.updZDot()[++ix] = t_G[1];
            state.updZDot()[++ix] = t_G[2];

            // tDot
            state.updZDot()[++ix] = kn * n_G[0];
            state.updZDot()[++ix] = kn * n_G[1];
            state.updZDot()[++ix] = kn * n_G[2];

            // nDot
            state.updZDot()[++ix] = -kn * t_G[0] - tau_g * b_G[0];
            state.updZDot()[++ix] = -kn * t_G[1] - tau_g * b_G[1];
            state.updZDot()[++ix] = -kn * t_G[2] - tau_g * b_G[2];

            // bDot
            state.updZDot()[++ix] = tau_g * n_G[0];
            state.updZDot()[++ix] = tau_g * n_G[1];
            state.updZDot()[++ix] = tau_g * n_G[2];

            return 0;
        }

    private:
        // Surface geometry for the geodesic.
        ContactGeometry m_geometry;
        // Surface to ground transform.
        Transform m_X_GS;
        // Initial Frenet frame of the geodesic.
        Transform m_X_GP;
    };
};

UnconstrainedGeodesic::UnconstrainedGeodesic(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    const Transform& X_GP)
{
    adoptSystemGuts(new Guts(geometry, X_GS, X_GP));
    DefaultSystemSubsystem defsub(*this);
}

/* Compute the UnconstrainedGeodesic that should match the curve segment on an
obstacle within a CableSpan. Returns the initial and final state of the
unconstrained geodesic. */
std::array<UnconstrainedGeodesicState, 2>
calcUnconstrainedGeodesicAtCableObstacle(
    const State& state,
    const CableSubsystem& subsystem,
    CableSpanIndex cableIx,
    CableSpanObstacleIndex obsIx,
    const IntegratorTestParameters& parameters)
{
    SimTK_ASSERT_ALWAYS(
        subsystem.getCable(cableIx).isInContactWithObstacle(state, obsIx),
        "Curve segment does not exist if obstacle not in contact with cable");
    const Real finalArcLength =
        subsystem.getCable(cableIx).calcCurveSegmentArcLength(state, obsIx);

    // Create the UnconstrainedGeodesic system.
    UnconstrainedGeodesic sys(state, subsystem, cableIx, obsIx);

    // Setup the integrator for computing the unconstrained geodesic.
    RungeKuttaMersonIntegrator integ(sys);
    integ.setAccuracy(parameters.accuracy);
    integ.setFinalTime(finalArcLength);

    // Initialize the unconstrained geodesic state.
    State initState = sys.realizeTopology();
    initState.setTime(0.);
    integ.initialize(initState);

    // Hold on to the initial state; We return it at the function end.
    UnconstrainedGeodesicState y0 =
        UnconstrainedGeodesic::getGeodesicState(integ.getState());

    // Integrate to fhe final arc length.
    constexpr int c_maxIntegrationSteps = 1000000;
    for (int i = 0; i < c_maxIntegrationSteps; ++i) {
        if (integ.stepTo(finalArcLength) == Integrator::EndOfSimulation) {
            break;
        }
    }
    SimTK_ASSERT_ALWAYS(
        std::abs(integ.getState().getTime() - finalArcLength) < Eps,
        "Failed to integrate UnconstrainedGeodesicState to final arc length");

    // Get the state at the final arc length.
    UnconstrainedGeodesicState y1 =
        UnconstrainedGeodesic::getGeodesicState(integ.getState());

    // Return the initial and final state of the UnconstrainedGeodesic.
    return {y0, y1};
}

} // namespace

//==============================================================================
//                      Geodesic test for the CableSubsystem
//==============================================================================
namespace
{

// Given a CableSpan's curve segment compute the unconstrained geodesic from
// the same initial conditions and assert correctness of both the unconstrained
// geodesic and the curve segment.
void runCurveSegmentGeodesicTest(
    const State& state,
    const CableSubsystem& subsystem,
    CableSpanIndex cableIx,
    CableSpanObstacleIndex obsIx,
    const IntegratorTestParameters& parameters,
    std::ostream& os)
{
    const CableSpan& cable = subsystem.getCable(cableIx);

    // Skip obstacles that are not in contact with the cable.
    if (!cable.isInContactWithObstacle(state, obsIx)) {
        os << "skip Geodesic Test of";
        os << ": cable " << cableIx;
        os << ", obstacle " << obsIx;
        os << ": Obstacle not in contact\n";
        return;
    }

    os << "START Geodesic Test of";
    os << ": cable " << cableIx;
    os << ", obstacle " << obsIx;
    os << ", with length " << cable.calcCurveSegmentArcLength(state, obsIx)
       << "\n";

    // First some sanity checks on the configuration parameters.

    // The accuracy of the unconstrained geodesic state should be relatively
    // small (= very accurate). How small depends on the model scale, but let's
    // just make sure it is not too large here.
    constexpr Real c_minAccuracy = 1e-10;
    SimTK_ERRCHK2_ALWAYS(
        parameters.accuracy <= c_minAccuracy,
        "CableSubsystemTestHelper::Impl::runGeodesicTest",
        "Invalid parameter: Unconstrained geodesic accuracy (=%e) should be "
        "less or equal to %e",
        parameters.accuracy,
        c_minAccuracy);
    SimTK_ERRCHK2_ALWAYS(
        parameters.accuracy < parameters.driftTolerance,
        "CableSubsystemTestHelper::Impl::runGeodesicTest",
        "Invalid parameter: Unconstrained geodesic accuracy (=%e) should be "
        "smaller than the drift tolerance (=%e)",
        parameters.accuracy,
        parameters.driftTolerance);
    SimTK_ERRCHK2_ALWAYS(
        parameters.accuracy < parameters.assertionTolerance,
        "CableSubsystemTestHelper::Impl::runGeodesicTest",
        "Invalid parameter: Unconstrained geodesic accuracy (=%e) should be "
        "smaller than the assertion tolerance (=%e)",
        cable.getCurveSegmentAccuracy(),
        parameters.assertionTolerance);
    SimTK_ERRCHK2_ALWAYS(
        cable.getCurveSegmentAccuracy() < parameters.assertionTolerance,
        "CableSubsystemTestHelper::Impl::runGeodesicTest",
        "Invalid parameter: Curve segment accuracy (=%e) should be smaller "
        "than the assertion tolerance (=%e)",
        cable.getCurveSegmentAccuracy(),
        parameters.assertionTolerance);

    // Compute the unconstrained geodesic that should match the cable obstacle's
    // curve segment.
    std::array<UnconstrainedGeodesicState, 2> unconstrainedGeodesicStates =
        calcUnconstrainedGeodesicAtCableObstacle(
            state,
            subsystem,
            cableIx,
            obsIx,
            parameters);
    // Alias for the initial state.
    const UnconstrainedGeodesicState& y0 = unconstrainedGeodesicStates.front();
    // Alias for the final state.
    const UnconstrainedGeodesicState& y1 = unconstrainedGeodesicStates.back();

    // Helper for computing the surface constraint value given the
    // UnconstrainedGeodesicState.
    auto calcSurfaceValue = [&](const UnconstrainedGeodesicState& y) -> Real
    {
        const ContactGeometry& geometry =
            subsystem.getCable(cableIx).getObstacleContactGeometry(obsIx);
        const Transform X_GS =
            getXformSurfaceToGround(state, subsystem, cableIx, obsIx);
        const Vec3 x_S = X_GS.shiftBaseStationToFrame(y.x);
        return geometry.calcSurfaceValue(x_S);
    };

    // Helper for computing the surface unit normal given the
    // UnconstrainedGeodesicState.
    auto calcSurfaceUnitNormal =
        [&](const UnconstrainedGeodesicState& y) -> UnitVec3
    {
        const ContactGeometry& geometry =
            subsystem.getCable(cableIx).getObstacleContactGeometry(obsIx);
        const Transform X_GS =
            getXformSurfaceToGround(state, subsystem, cableIx, obsIx);
        const Vec3 x_S = X_GS.shiftBaseStationToFrame(y.x);
        return UnitVec3(
            X_GS.xformFrameVecToBase(geometry.calcSurfaceUnitNormal(x_S)));
    };

    // This max error we only use to print some information for the user.
    Real maxErr = 0.;
    // Helper for updating the max error, given a new error evaluation.
    auto trackMaxErr = [&](Real err) -> Real
    {
        maxErr = std::max(err, maxErr);
        return err;
    };

    // Assert the initial state of the UnconstrainedGeodesic.
    SimTK_ASSERT_ALWAYS(
        trackMaxErr(std::abs(calcSurfaceValue(y0))) <
            parameters.constraintTolerance,
        "Initial point of UnconstrainedGeodesic not on surface");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y0.R.col(0).norm() - 1.)) < parameters.constraintTolerance,
        "Initial tangent of UnconstrainedGeodesic does not have unit norm");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y0.R.col(1).norm() - 1.)) < parameters.constraintTolerance,
        "Initial normal of UnconstrainedGeodesic does not have unit norm");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y0.R.col(2).norm() - 1.)) < parameters.constraintTolerance,
        "Initial binormal of UnconstrainedGeodesic does not have unit norm");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((calcSurfaceUnitNormal(y0) - y0.R.col(1)).norm()) <
            parameters.constraintTolerance,
        "Initial normal of UnconstrainedGeodesic does not match surface normal");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y0.R * (~y0.R) - Mat33(1)).norm()) <
            parameters.constraintTolerance,
        "Initial axes of UnconstrainedGeodesic are not mutually perpendicular");
    os << "PASSED TEST: Initial UnconstrainedGeodesicState on surface";
    os << " (err = " << maxErr << " < " << parameters.constraintTolerance
       << " = eps)" << std::endl;
    maxErr = 0.; // Reset the max error after each section of tests.

    // Assert the final state of the UnconstrainedGeodesic.
    SimTK_ASSERT_ALWAYS(
        trackMaxErr(std::abs(calcSurfaceValue(y1))) < parameters.driftTolerance,
        "Final point of UnconstrainedGeodesic not on surface");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y1.R.col(0).norm() - 1.)) < parameters.driftTolerance,
        "Final tangent of UnconstrainedGeodesic does not have unit norm");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y1.R.col(1).norm() - 1.)) < parameters.driftTolerance,
        "Final normal of UnconstrainedGeodesic does not have unit norm");
    SimTK_ASSERT1_ALWAYS(
        trackMaxErr((y1.R.col(2).norm() - 1.)) < parameters.driftTolerance,
        "Final binormal of UnconstrainedGeodesic does not have unit norm (%e)",
        y1.R.col(2).norm() - 1.);
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((calcSurfaceUnitNormal(y1) - y1.R.col(1)).norm()) <
            parameters.driftTolerance,
        "Final normal of UnconstrainedGeodesic does not match surface normal");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((y1.R * (~y1.R) - Mat33(1)).norm()) <
            parameters.driftTolerance,
        "Final axes of UnconstrainedGeodesic are not mutually perpendicular");
    os << "PASSED TEST: Final UnconstrainedGeodesic on surface";
    os << " (err = " << maxErr << " < " << parameters.driftTolerance
       << " = eps)" << std::endl;
    maxErr = 0.;

    // Assert that the curve segment matches the UnconstrainedGeodesic at the
    // first contact point.
    const Transform X_GP =
        cable.calcCurveSegmentInitialFrenetFrame(state, obsIx);
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((X_GP.p() - y0.x).norm()) < parameters.assertionTolerance,
        "Curve segment Frenet frame at initial contact point does not match UnconstrainedGeodesic");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((X_GP.R().asMat33() * (~y0.R) - Mat33(1)).norm()) <
            parameters.assertionTolerance,
        "Curve segment Frenet frame at initial contact point does not match UnconstrainedGeodesic");
    os << "PASSED TEST: Initial curve segment frame on surface";
    os << " (err = " << maxErr << " < " << parameters.assertionTolerance
       << " = eps)" << std::endl;
    maxErr = 0.;

    // Assert that the curve segment matches the UnconstrainedGeodesic at the
    // final contact point.
    const Transform X_GQ = cable.calcCurveSegmentFinalFrenetFrame(state, obsIx);
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((X_GQ.p() - y1.x).norm()) < parameters.assertionTolerance,
        "Curve segment Frenet frame at final contact point does not match UnconstrainedGeodesic");
    SimTK_ASSERT_ALWAYS(
        trackMaxErr((X_GQ.R().asMat33() * (~y1.R) - Mat33(1)).norm()) <
            parameters.assertionTolerance,
        "Curve segment Frenet frame at final contact point does not match UnconstrainedGeodesic");
    os << "PASSED TEST: Final curve segment frame matches "
          "UnconstrainedGeodesic";
    os << " (err = " << maxErr << " < " << parameters.assertionTolerance
       << " = eps)" << std::endl;
}

} // namespace

void CableSubsystemTestHelper::Impl::runGeodesicTest(
    const State& state,
    const CableSubsystem& subsystem,
    std::ostream& testReport) const
{
    // Call the geodesic test for each curve segment of each cable in the
    // subsystem.
    for (CableSpanIndex cableIx(0); cableIx < subsystem.getNumCables();
         ++cableIx) {
        for (CableSpanObstacleIndex obsIx(0);
             obsIx < subsystem.getCable(cableIx).getNumObstacles();
             ++obsIx) {
            runCurveSegmentGeodesicTest(
                state,
                subsystem,
                cableIx,
                obsIx,
                m_integratorTestParameters,
                testReport);
        }
    }
}

//==============================================================================
//                  Implementation of CableSubsystemTestHelper
//==============================================================================

void CableSubsystemTestHelper::Impl::testCurrentPath(
    const State& state,
    const CableSubsystem& subsystem,
    std::ostream& testReport) const
{
    // Intentionally skip the perturbation test if the geodesic test fails.
    runGeodesicTest(state, subsystem, testReport);
    runPerturbationTest(state, subsystem, testReport);
}

CableSubsystemTestHelper::CableSubsystemTestHelper() :
    m_impl(new CableSubsystemTestHelper::Impl())
{}

CableSubsystemTestHelper::~CableSubsystemTestHelper() noexcept = default;

CableSubsystemTestHelper::CableSubsystemTestHelper(
    const CableSubsystemTestHelper& source) :
    m_impl(new CableSubsystemTestHelper::Impl(*source.m_impl))
{}

CableSubsystemTestHelper& CableSubsystemTestHelper::operator=(
    const CableSubsystemTestHelper& source)
{
    return *this = CableSubsystemTestHelper(source);
}

CableSubsystemTestHelper::CableSubsystemTestHelper(
    CableSubsystemTestHelper&& source) noexcept = default;

CableSubsystemTestHelper& CableSubsystemTestHelper::operator=(
    CableSubsystemTestHelper&& source) noexcept = default;

void CableSubsystemTestHelper::testCurrentPath(
    const State& state,
    const CableSubsystem& subsystem,
    std::ostream& testReport) const
{
    m_impl->testCurrentPath(state, subsystem, testReport);
}
