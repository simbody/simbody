/* -------------------------------------------------------------------------- *
 *            Simbody(tm) Example: Test Cable Over Smooth Surfaces            *
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

// This example outlines a unit test for a CableSpan. The path error jacobian and
// the individual geodesics are tested for a cable spanning two surfaces.

#include "SimTKcommon/internal/Transform.h"
#include "Simbody.h"
#include "simbody/internal/CableSpan.h"
#include "simmath/internal/ContactGeometry.h"
#include <cassert>
#include <iostream>
#include <ostream>
#include <stdexcept>

using namespace SimTK;
using FrenetFrame = Transform;
static const CoordinateAxis TangentAxis  = CoordinateAxis::XCoordinateAxis();
static const CoordinateAxis NormalAxis   = CoordinateAxis::YCoordinateAxis();
static const CoordinateAxis BinormalAxis = CoordinateAxis::ZCoordinateAxis();

//==============================================================================
//                      GEODESIC FORWARD INTEGRATION
//==============================================================================

// TODO Sign of implicit function was flipped.
Vec3 calcSurfaceConstraintGradient(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    Vec3 point_G)
{
    std::cout << "WARNING: using local calcSurfaceConstraintGradient function\n";
    std::cout << "WARNING: FLIPPED SIGN\n";
    const Vec3 point_S = X_GS.shiftBaseStationToFrame(point_G);
    return -X_GS.xformFrameVecToBase(geometry.calcSurfaceGradient(point_S));
}

// TODO Sign of implicit function was flipped.
Mat33 calcSurfaceConstraintHessian(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    Vec3 point_G)
{
    std::cout << "WARNING: using local calcSurfaceConstraintHessian function\n";
    std::cout << "WARNING: FLIPPED SIGN\n";
    const Vec3 point_S = X_GS.shiftBaseStationToFrame(point_G);
    return -X_GS.R() * geometry.calcSurfaceHessian(point_S) * X_GS.RInv();
}

Real calcNormalCurvature(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    Vec3 point,
    Vec3 tangent)
{
    // TODO use:
    // const UnitVec3 n = geometry.calcSurfaceUnitNormal(x_S);
    std::cout << "WARNING: using local calcNormalCurvature function\n";
    const Vec3& p  = point;
    const Vec3& v  = tangent;
    const Vec3 g   = calcSurfaceConstraintGradient(geometry, X_GS, p);
    const Vec3 h_v = calcSurfaceConstraintHessian(geometry, X_GS, p) * v;
    // Sign flipped compared to thesis: kn = negative, see eq 3.63
    return -dot(v, h_v) / g.norm();
}

UnitVec3 calcSurfaceNormal(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    Vec3 point_G)
{
    std::cout << "WARNING: using local calcSurfaceNormal function\n";
    // TODO use:
    // const Real kn = geometry.calcSurfaceCurvatureInDirection(x_S, t_S);
    return UnitVec3{calcSurfaceConstraintGradient(geometry, X_GS, point_G)};
}

Real calcGeodesicTorsion(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    Vec3 point_G,
    Vec3 tangent_G)
{
    // TODO dont do this.
    std::cout << "WARNING: using local calcGeodesicTorsion function\n";
    const Vec3& p  = point_G;
    const Vec3& v  = tangent_G;
    const Vec3 g   = calcSurfaceConstraintGradient(geometry, X_GS, point_G);
    const Vec3 h_v = calcSurfaceConstraintHessian(geometry, X_GS, point_G) * v;
    const Vec3 gxv = cross(g, v);
    return -dot(h_v, gxv) / dot(g, g);
}

struct GeodesicState
{
    static constexpr int N = 16;

    using Derivative = Vec<N, Real>;

    const Vec<N, Real>& asVec() const
    {
        return reinterpret_cast<const Vec<N, Real>&>(xDot[0]);
    }

    Vec<N, Real>& asVecMut()
    {
        return reinterpret_cast<Vec<N, Real>&>(xDot[0]);
    }

    static Derivative CalcDerivative(
        const ContactGeometry& geometry,
        const Transform& X_GS,
        const GeodesicState& q)
    {
        GeodesicState qDot;

        const Real kn    = calcNormalCurvature(geometry, X_GS, q.x, q.t);
        const UnitVec3 n = calcSurfaceNormal(geometry, X_GS, q.x);
        const Real tau_g = calcGeodesicTorsion(geometry, X_GS, q.x, q.t);

        qDot.xDot = kn * n;
        qDot.x    = q.xDot;

        qDot.t = kn * n;
        qDot.n = -kn * q.t + tau_g * q.b;
        qDot.b = -tau_g * q.n;

        qDot.l = 1.;

        return qDot.asVec();
    }

    Vec3 xDot = {NAN, NAN, NAN};
    Vec3 x    = {NAN, NAN, NAN};
    Vec3 t    = {NAN, NAN, NAN};
    Vec3 n    = {NAN, NAN, NAN};
    Vec3 b    = {NAN, NAN, NAN};
    Real l    = 0.;
};

GeodesicState operator-(const GeodesicState& lhs, const GeodesicState& rhs)
{
    GeodesicState out;
    out.asVecMut() = lhs.asVec() - rhs.asVec();
    return out;
}

GeodesicState operator+(
    const GeodesicState& lhs,
    const Vec<GeodesicState::N, Real>& rhs)
{
    GeodesicState out;
    out.asVecMut() = lhs.asVec() + rhs;
    return out;
}

//==============================================================================
//                      RUNGE KUTTA 4
//==============================================================================

void RungeKutta4Step(
    const ContactGeometry geometry,
    const Transform& X_GS,
    GeodesicState& y,
    Real dl)
{
    using Y = GeodesicState;
    using D = GeodesicState::Derivative;

    auto f = [&](const Y& yk) -> D {
        return Y::CalcDerivative(geometry, X_GS, yk);
    };

    D k0, k1, k2, k3;

    {
        const Y& yk = y;
        k0          = f(yk);
    }

    {
        const Real h = dl / 2.;
        Y yk         = y + (h * k0);
        k1           = f(yk);
    }

    {
        const Real h = dl / 2.;
        Y yk         = y + (h * k1);
        k2           = f(yk);
    }

    {
        const Real h = dl;
        Y yk         = y + (h * k2);
        k3           = f(yk);
    }

    const Real w = dl / 6.;

    y = y + (w * k0 + (w * 2.) * k1 + (w * 2.) * k2 + w * k3);
}

std::vector<GeodesicState> RecomputeUnconstrainedGeodesic(
    const State& state,
    const CurveSegment& curve,
    size_t integratorSteps)
{
    std::vector<GeodesicState> samples;

    const FrenetFrame& K_P = curve.getFrenetFrameStart(state);

    GeodesicState q;

    q.xDot = Vec3(K_P.R().getAxisUnitVec(TangentAxis));
    q.x    = K_P.p();
    q.t    = Vec3(K_P.R().getAxisUnitVec(TangentAxis));
    q.n    = Vec3(K_P.R().getAxisUnitVec(NormalAxis));
    q.b    = Vec3(K_P.R().getAxisUnitVec(BinormalAxis));
    q.l    = 0.;

    samples.push_back(q);

    const int n              = integratorSteps;
    const ContactGeometry& g = curve.getContactGeometry();
    const Mobod& body        = curve.getMobilizedBody();
    const Transform X_GS =
        body.getBodyTransform(state).compose(curve.getXformSurfaceToBody());

    Real l  = curve.getSegmentLength(state);
    Real dl = l / static_cast<Real>(n);

    for (int k = 0; k < n; ++k) {
        RungeKutta4Step(g, X_GS, q, dl);
        samples.push_back(q);
    }

    if (std::abs(l - samples.back().l) > Eps * 1e2) {
        throw std::runtime_error(
            "Incorrect geodesic length after integration.");
    }

    return samples;
}

//==============================================================================
//                      UNIT TESTING
//==============================================================================

bool AssertRecomputedPosition(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    const GeodesicState& q,
    const std::string& msg,
    std::ostream& os,
    Real eps)
{
    const Vec3 x_S = X_GS.shiftBaseStationToFrame(q.x);
    const Real c   = geometry.calcSurfaceValue(x_S);

    if (std::abs(c) > eps) {
        os << "TEST FAILED at " << msg << " (l = " << q.l << ")\n";

        os << "--> Cause:             ";
        os << "Recomputed position of " << msg << " has left the surface\n";

        os << "--> Failed assertion:  ";
        os << "|ContactGeometry::calcSurfaceValue(" << x_S << ")|";
        os << " = |" << c << "|";
        os << " < " << eps << "\n";

        os << "--> Possible errors in:\n";
        os << "1. The surface to ground transform.\n";
        os << "2. The surface constraint value.\n";
        os << "3. The normal curvature.\n";
        os << "4. The surface normal vector.\n";

        return false;
    }

    return true;
}

bool AssertRecomputedFrenetFrame(
    const ContactGeometry& geometry,
    const Transform& X_GS,
    const GeodesicState& q,
    const std::string& msg,
    std::ostream& os,
    Real eps)
{
    std::ostringstream oss;

    auto AssertAxisNorm = [&](const Vec3& v,
                              const std::string& whichAxis) -> bool {
        if (std::abs(v.norm() - 1.) > eps) {
            oss << "--> Cause: ";
            oss << whichAxis << " does not have unit norm";
            oss << " (eps = " << eps << ")\n";

            oss << "--> Failed assertion: ";
            oss << "|" << whichAxis << "|";
            oss << " = |" << v << "|";
            oss << " = " << v.norm();
            oss << " = 1 \n";

            return false;
        }
        return true;
    };

    auto AssertOrthogonality = [&](const Vec3& a,
                                   const Vec3& b,
                                   const Vec3& c,
                                   const std::string& whichAxisIsA,
                                   const std::string& whichAxisIsB,
                                   const std::string& whichAxisIsC) -> bool {
        bool success = true;

        const Real dotErr = std::abs(dot(a, b));
        if (dotErr > eps) {
            oss << "--> Cause: ";
            oss << whichAxisIsA << " and " << whichAxisIsB
                << " are not perpendicular";
            oss << " (eps = " << eps << ")\n";

            oss << "--> Failed assertion: ";
            oss << "|dot(" << whichAxisIsA << ", " << whichAxisIsB << ")|";
            oss << " = |dot(" << a << ", " << b << ")|";
            oss << " = " << dotErr << " < " << eps << "\n";

            success = false;
        }

        const Real crossErr = (cross(a, b) - c).norm();
        if (crossErr > eps) {
            oss << "--> Cause: ";
            oss << "Cross product of " << whichAxisIsA << " and "
                << whichAxisIsB << " should give " << whichAxisIsC << "\n";

            oss << "--> Failed assertion: ";
            oss << "|" << whichAxisIsA << " % " << whichAxisIsB << " - "
                << whichAxisIsC << "| = ";
            oss << "|" << a << " % " << b << " - ( " << c
                << " ) | = " << crossErr << " < " << eps << "\n";
            success = false;
        }

        return success;
    };

    bool success = true;

    success &= AssertAxisNorm(q.t, "Tangent");
    success &= AssertAxisNorm(q.n, "Normal");
    success &= AssertAxisNorm(q.b, "Binormal");

    success &=
        AssertOrthogonality(q.t, q.n, q.b, "Tangent", "Normal", "Binormal");
    success &=
        AssertOrthogonality(q.b, q.t, q.n, "Tangent", "Binormal", "Normal");
    success &=
        AssertOrthogonality(q.n, q.b, q.t, "Normal", "Binormal", "Tangent");

    if ((q.xDot - q.t).norm() > eps) {
        oss << "--> Cause: ";
        oss << "Position derivative should be equal to tangent\n";

        oss << "--> Failed assertion:";
        oss << "|xDot - t|";
        oss << " = |" << q.xDot << " - ( " << q.t << " )|";
        oss << " = " << (q.xDot - q.t);
        oss << " < " << eps << "\n";
    }

    if (!success) {
        os << "TEST FAILED at " << msg << " (l = " << q.l << ")\n";
        os << oss.str();
        os << "--> Possible errors in computing the geodesic torsion.\n";
    }

    return success;
}

bool AssertCurveSegmentFrameEqualsRecomputedFrame(
    const FrenetFrame& K,
    const GeodesicState& q,
    const std::string& msg,
    std::ostream& os,
    Real eps)
{
    std::ostringstream oss;

    auto AssertEqual =
        [&](const Vec3& a, const Vec3& b, const std::string& desc) -> bool {
        if ((a - b).norm() > eps) {
            oss << "--> Cause: ";
            oss << desc << " of curve segment does not match recomputed "
                << desc;
            oss << " (at length = " << q.l << ")\n";

            oss << "--> Failed assertion: ";
            oss << "|" << a << " - (" << b << ")|";
            oss << " = " << (a - b).norm();
            oss << " < " << eps << "\n";

            return false;
        }
        return true;
    };

    bool success = true;

    success &= AssertEqual(K.p(), q.x, "Position");
    success &= AssertEqual(K.R().getAxisUnitVec(TangentAxis), q.t, "Tangent");
    success &= AssertEqual(K.R().getAxisUnitVec(NormalAxis), q.n, "Normal");
    success &= AssertEqual(K.R().getAxisUnitVec(BinormalAxis), q.b, "Binormal");

    if (!success) {
        os << "TEST FAILED at " << msg << "\n";
        os << oss.str();
        os << "--> Possible errors in the geodesic computation within the "
              "CurveSegment\n";
    }

    return success;
}

bool RunRungeKutta4Test(std::ostream& os)
{
    os << "Start all wrapping tests";

    return true;
}

bool RunCurveSegmentShooterTest(
    const State& state,
    const CurveSegment& curve,
    Real eps,
    int integratorSteps,
    std::ostream& os)
{

    std::vector<GeodesicState> testSamples =
        RecomputeUnconstrainedGeodesic(state, curve, integratorSteps);

    const ContactGeometry& g = curve.getContactGeometry();
    const Mobod& body        = curve.getMobilizedBody();
    const Transform X_GS =
        body.getBodyTransform(state).compose(curve.getXformSurfaceToBody());

    bool success = true;
    success &= AssertRecomputedPosition(
        g,
        X_GS,
        testSamples.front(),
        "first contact point",
        os,
        eps);
    success &= AssertRecomputedPosition(
        g,
        X_GS,
        testSamples.back(),
        "last contact point",
        os,
        eps);

    for (const GeodesicState& q : testSamples) {
        success =
            success &&
            AssertRecomputedPosition(g, X_GS, q, "interpolated point", os, eps);
    }

    for (const GeodesicState& q : testSamples) {
        success = success && AssertRecomputedFrenetFrame(
                                 g,
                                 X_GS,
                                 q,
                                 "unconstrained frenet frame.",
                                 os,
                                 eps);
    }

    std::vector<Transform> frames;
    std::function<void(Vec3 point_G, UnitVec3 tangent_G)> Logger =
        [&](Vec3 point_G, UnitVec3 tangent_G) {
            Transform frame;
            frame.updP() = point_G;
            frame.updR().setRotationFromTwoAxes(
                calcSurfaceNormal(curve.getContactGeometry(), X_GS, point_G),
                NormalAxis,
                tangent_G,
                TangentAxis);
            frames.push_back(frame);
        };

    curve.calcPathPointsAndTangents(state, Logger, integratorSteps + 1);

    if (frames.size() != testSamples.size()) {
        std::cout << frames.size() << "\n";
        std::cout << testSamples.size() << "\n";
        throw std::runtime_error("Incorrect number of samples");
    }

    success = success && AssertCurveSegmentFrameEqualsRecomputedFrame(
                             frames.front(),
                             testSamples.front(),
                             "Equals unconstrained at start",
                             os,
                             eps);

    success = success && AssertCurveSegmentFrameEqualsRecomputedFrame(
                             frames.back(),
                             testSamples.back(),
                             "Equals unconstrained at end",
                             os,
                             eps);

    for (size_t i = 0; i < frames.size(); ++i) {
        success = success && AssertCurveSegmentFrameEqualsRecomputedFrame(
                                 frames.at(i),
                                 testSamples.at(i),
                                 "Equals unconstrained at interpolated",
                                 os,
                                 eps);
    }

    return success;
}

int main()
{
    const Real integratorAccuracy = 1e-6;
    const Real projectionAccuracy = 1e-6;
    const int projectionMaxIter   = 100;
    const int pathMaxIter         = 100;
    const Real pathAccuracy       = 1e-4;
    const Real maxStepDeg         = 5.;

    try {
        // Create the system.
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);
        GeneralForceSubsystem forces(system);

        system.setUseUniformBackground(true); // no ground plane in display

        Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));

        Body::Rigid ballBody(MassProperties(1.0, Vec3(0), Inertia(1)));
        const Real Rad = 1.;

        MobilizedBody Ground = matter.Ground();

        MobilizedBody::Ball ball0(
            Ground,
            Transform(Vec3(0)),
            ballBody,
            Transform(Vec3(0)));

        MobilizedBody::Ball ball1(
            Ground,
            Transform(Vec3(0)),
            ballBody,
            Transform(Vec3(0)));

        CableSpan cable(
            cables,
            Ground,
            Vec3(-2., 0.1, 0.), // origin
            Ground,
            Vec3(1.75, 0.05, 0.)); // termination

        cable.setSurfaceConstraintTolerance(projectionAccuracy);
        cable.setSurfaceProjectionMaxIter(projectionMaxIter);
        cable.setPathErrorAccuracy(pathAccuracy);
        cable.setSolverMaxIter(pathMaxIter);
        cable.setIntegratorAccuracy(integratorAccuracy);
        cable.setMaxRadialStepInDegrees(maxStepDeg);

        cable.addSurfaceObstacle(
            ball0,
            Transform(Vec3{-0.5, 0., 0.}),
            ContactGeometry::Sphere(Rad),
            {0., 1., 0.});

        cable.addSurfaceObstacle(
            ball1,
            Transform(Vec3{1.1, 0.2, 0.3}),
            ContactGeometry::Ellipsoid(Vec3{0.5, 0.7, 1.2} * Rad),
            {0., 1., 0.});

        Visualizer viz(system);
		viz.setShowFrameNumber(true);
		system.addEventReporter(new Visualizer::Reporter(viz, 0.1 * 1. / 30));

        // Initialize the system and state.
        system.realizeTopology();
        State s = system.getDefaultState();

        system.realize(s, Stage::Position);
		viz.report(s);

        std::ostringstream oss;
        bool testPassed = true;
        for (CurveSegmentIndex ix = CurveSegmentIndex(0);
             ix < cable.getNumSurfaceObstacles();
             ++ix) {
            testPassed &= RunCurveSegmentShooterTest(
                s,
                cable.getCurveSegment(ix),
                cable.getIntegratorAccuracy(),
                1000,
                oss);
        }

        CableSubsystemTestHelper perturbationTestHelper;
        testPassed &= perturbationTestHelper
                          .applyPerturbationTest(s, cables, 1e-4, 1e-3, oss);

        std::cout << oss.str() << "\n";

        if (!testPassed) {
            throw std::runtime_error("Failed test");
        }

        std::cout << "All tests passed\n";
        // TODO add more configurations etc..

    } catch (const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << "\n";
        return -1.;
    }
    return 0;
}
