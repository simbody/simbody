/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-26 Stanford University and the Authors.        *
 * Authors: Nicholas Bianco                                                   *
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

#include "SimTKsimbody.h"

using namespace SimTK;

namespace {

    // A branched multibody of bodies containing various mobilizers of different
    // types including Pin, Ellipsoid, and CantileverFreeBeam mobilizers.
    class BranchedSystem {
    public:
        enum MobilizerType {Ground=0, Pin=1, Ellipsoid=2, CantileverFreeBeam=3,
                            Gimbal=4, Weld=5, Free=6};

        BranchedSystem() :
                m_matter(m_system), m_forces(m_system),
                m_gravity(m_forces, m_matter, -YAxis, 9.8) {

            Body::Rigid body(MassProperties(1.0, Vec3(0), UnitInertia(1)));

            m_pin = MobilizedBody::Pin(
                    m_matter.Ground(), X_PF[Pin],
                    body, X_BM[Pin]);

            // Branch 1
            m_ellipsoid = MobilizedBody::Ellipsoid(
                    m_pin, X_PF[Ellipsoid],
                    body, X_BM[Ellipsoid],
                    m_ellipsoidRadii);

            m_cantileverFreeBeam = MobilizedBody::CantileverFreeBeam(
                    m_ellipsoid, X_PF[CantileverFreeBeam],
                    body, X_BM[CantileverFreeBeam],
                    m_cantileverFreeBeamLength);

            // Branch 2
            m_gimbal = MobilizedBody::Gimbal(
                    m_pin, X_PF[Gimbal],
                    body, X_BM[Gimbal]);

            m_weld = MobilizedBody::Weld(
                    m_gimbal, X_PF[Weld],
                    body, X_BM[Weld]);

            m_free = MobilizedBody::Free(
                    m_weld, X_PF[Free],
                    body, X_BM[Free]);
        }

        void loadDefaultState(State& state) {
            for (int i = 0; i < state.getNQ(); ++i) {
                state.updQ()[i] = 0.1 * (i+1);
            }
            for (int i = 0; i < state.getNU(); ++i) {
                state.updU()[i] = 1.0 * (i+1);
            }
        }

        MultibodySystem        m_system;
        SimbodyMatterSubsystem m_matter;
        GeneralForceSubsystem  m_forces;
        Force::Gravity         m_gravity;

        MobilizedBody::Pin                m_pin;
        MobilizedBody::Ellipsoid          m_ellipsoid;
        MobilizedBody::CantileverFreeBeam m_cantileverFreeBeam;
        MobilizedBody::Gimbal             m_gimbal;
        MobilizedBody::Weld               m_weld;
        MobilizedBody::Free               m_free;

        const Real m_cantileverFreeBeamLength = 1.23;
        const Vec3 m_ellipsoidRadii = Vec3(0.1, 0.2, 0.3);

    private:
        const Array_<Transform> X_PF = {
            Transform(),                              // [0] Ground
            Transform(Rotation(BodyRotationSequence,  // [1] Pin
                            Pi/8, XAxis,
                            Pi/7, YAxis),
                    Vec3(-0.4, 0.5, -0.6)),
            Transform(Rotation(BodyRotationSequence,  // [2] Ellipsoid
                            Pi/4, ZAxis,
                            Pi/5, XAxis),
                    Vec3(0.1, -0.2, 0.3)),
            Transform(Rotation(BodyRotationSequence,  // [3] CantileverFreeBeam
                            -Pi/3, XAxis,
                            Pi/5, ZAxis),
                    Vec3(0.6, -0.7, 0.8)),
            Transform(Rotation(BodyRotationSequence,  // [4] Gimbal
                            Pi/6, YAxis,
                            -Pi/7, ZAxis),
                    Vec3(-0.9, 1.0, -1.1)),
            Transform(Rotation(BodyRotationSequence,  // [5] Weld
                            Pi/8, ZAxis,
                            Pi/9, YAxis),
                    Vec3(1.2, -1.3, 1.4)),
            Transform(Rotation(BodyRotationSequence,  // [6] Free
                            -Pi/10, YAxis,
                            Pi/11, XAxis),
                    Vec3(-1.5, 1.6, -1.7))};
        const Array_<Transform> X_BM = {
            Transform(),                              // [0] Ground
            Transform(Rotation(BodyRotationSequence,  // [1] Pin
                            Pi/5, ZAxis,
                            Pi/6, YAxis),
                    Vec3(-0.10, 0.11, -0.12)),
            Transform(Rotation(BodyRotationSequence,  // [2] Ellipsoid
                            -Pi/4, XAxis,
                            Pi/3, YAxis),
                    Vec3(0.7, -0.8, 0.9)),
            Transform(Rotation(BodyRotationSequence,  // [3] CantileverFreeBeam
                            -Pi/6, ZAxis,
                            Pi/4, XAxis),
                    Vec3(-0.16, 0.17, -0.18)),
            Transform(Rotation(BodyRotationSequence,  // [4] Gimbal
                            Pi/5, YAxis,
                            -Pi/6, ZAxis),
                    Vec3(0.19, -0.20, 0.21)),
            Transform(Rotation(BodyRotationSequence,  // [5] Weld
                            -Pi/7, XAxis,
                            Pi/8, ZAxis),
                    Vec3(-0.22, 0.23, -0.24)),
            Transform(Rotation(BodyRotationSequence,  // [6] Free
                            Pi/9, ZAxis,
                            -Pi/10, YAxis),
                    Vec3(0.25, -0.26, 0.27))};
    };

}

void testMultiplyByPositionJacobianWrtInboardFramePositions() {
    BranchedSystem sys;
    State state = sys.m_system.realizeTopology();
    sys.loadDefaultState(state);
    sys.m_system.realize(state, Stage::Position);

    const int nb = sys.m_matter.getNumBodies();
    const Real h = 1e-5;

    Vector_<Vec3> dp_PF(nb);
    for (int b = 0; b < nb; ++b) {
        dp_PF[b] = Vec3(0.5*(b+1), -0.3*(b+1), 0.7*(b+1));
    }

    Vector_<Vec3> dp_GB;
    sys.m_matter.multiplyByPositionJacobianWrtInboardFramePositions(
            state, dp_PF, dp_GB);

    State pert = state;
    MobilizedBody mobods[6] = { sys.m_pin,
                                sys.m_ellipsoid,
                                sys.m_cantileverFreeBeam,
                                sys.m_gimbal,
                                sys.m_weld,
                                sys.m_free };
    for (int m = 0; m < 6; ++m) {
        const MobilizedBodyIndex bIdx = mobods[m].getMobilizedBodyIndex();
        Transform X_PF = mobods[m].getInboardFrame(pert);
        X_PF.updP() += h * dp_PF[bIdx];
        mobods[m].setInboardFrame(pert, X_PF);
    }
    sys.m_system.realize(pert, Stage::Position);

    for (int ib = 0; ib < nb; ++ib) {
        const Vec3 p0 = sys.m_matter.getMobilizedBody(MobilizedBodyIndex(ib))
                                    .getBodyTransform(state).p();
        const Vec3 p1 = sys.m_matter.getMobilizedBody(MobilizedBodyIndex(ib))
                                    .getBodyTransform(pert).p();
        SimTK_TEST_EQ_TOL(dp_GB[ib], (p1 - p0) / h, 1e-10);
    }
}

void testMultiplyByPositionJacobianWrtInboardFramePositionsTranspose() {
    BranchedSystem sys;
    State state = sys.m_system.realizeTopology();
    sys.loadDefaultState(state);
    sys.m_system.realize(state, Stage::Position);

    const int nb = sys.m_matter.getNumBodies();

    Vector_<Vec3> dp_PF(nb), g_GB(nb);
    for (int b = 0; b < nb; ++b) {
        dp_PF[b] = Vec3( 0.1*(b+1), -0.2*(b+1),  0.3*(b+1));
        g_GB[b]  = Vec3(-0.5*(b+1),  0.7*(b+1), -0.9*(b+1));
    }
    dp_PF[0] = Vec3(0);
    g_GB[0]  = Vec3(0);

    Vector_<Vec3> dp_GB, g_PF;
    sys.m_matter.multiplyByPositionJacobianWrtInboardFramePositions(
            state, dp_PF, dp_GB);
    sys.m_matter.multiplyByPositionJacobianWrtInboardFramePositionsTranspose(
            state, g_GB, g_PF);

    // <JpF*dp_PF, g_GB> = <dp_PF, ~JpF*g_GB>
    Real lhs = 0, rhs = 0;
    for (int b = 0; b < nb; ++b) {
        lhs += dot(g_GB[b],  dp_GB[b]);
        rhs += dot(dp_PF[b], g_PF[b]);
    }
    SimTK_TEST_EQ_TOL(lhs, rhs, 1e-10);
}

void testMultiplyByPositionJacobianWrtOutboardFramePositions() {
    BranchedSystem sys;
    State state = sys.m_system.realizeTopology();
    sys.loadDefaultState(state);
    sys.m_system.realize(state, Stage::Position);

    const int nb = sys.m_matter.getNumBodies();
    const Real h = 1e-5;

    Vector_<Vec3> dp_BM(nb);
    for (int b = 0; b < nb; ++b) {
        dp_BM[b] = Vec3(-0.2*(b+1), 0.9*(b+1), -0.4*(b+1));
    }

    Vector_<Vec3> dp_GB;
    sys.m_matter.multiplyByPositionJacobianWrtOutboardFramePositions(
        state, dp_BM, dp_GB);

    State pert = state;
    MobilizedBody mobods[6] = { sys.m_pin,
                                sys.m_ellipsoid,
                                sys.m_cantileverFreeBeam,
                                sys.m_gimbal,
                                sys.m_weld,
                                sys.m_free };
    for (int m = 0; m < 6; ++m) {
        const MobilizedBodyIndex bIdx = mobods[m].getMobilizedBodyIndex();
        Transform X_BM = mobods[m].getOutboardFrame(pert);
        X_BM.updP() += h * dp_BM[bIdx];
        mobods[m].setOutboardFrame(pert, X_BM);
    }
    sys.m_system.realize(pert, Stage::Position);

    for (int ib = 0; ib < nb; ++ib) {
        const Vec3 p0 = sys.m_matter.getMobilizedBody(MobilizedBodyIndex(ib))
                                    .getBodyTransform(state).p();
        const Vec3 p1 = sys.m_matter.getMobilizedBody(MobilizedBodyIndex(ib))
                                    .getBodyTransform(pert).p();
        SimTK_TEST_EQ_TOL(dp_GB[ib], (p1 - p0) / h, 1e-10);
    }
}

void testMultiplyByPositionJacobianWrtOutboardFramePositionsTranspose() {
    BranchedSystem sys;
    State state = sys.m_system.realizeTopology();
    sys.loadDefaultState(state);
    sys.m_system.realize(state, Stage::Position);

    const int nb = sys.m_matter.getNumBodies();

    Vector_<Vec3> dp_BM(nb), g_GB(nb);
    for (int b = 0; b < nb; ++b) {
        dp_BM[b] = Vec3( 0.2*(b+1), -0.4*(b+1),  0.6*(b+1));
        g_GB[b]  = Vec3(-0.3*(b+1),  0.9*(b+1), -0.5*(b+1));
    }
    dp_BM[0] = Vec3(0);
    g_GB[0]  = Vec3(0);

    Vector_<Vec3> dp_GB, g_BM;
    sys.m_matter.multiplyByPositionJacobianWrtOutboardFramePositions(
            state, dp_BM, dp_GB);
    sys.m_matter.multiplyByPositionJacobianWrtOutboardFramePositionsTranspose(
            state, g_GB, g_BM);

    // <JpM*dp_BM, g_GB> = <dp_BM, ~JpM*g_GB>
    Real lhs = 0, rhs = 0;
    for (int b = 0; b < nb; ++b) {
        lhs += dot(g_GB[b],  dp_GB[b]);
        rhs += dot(dp_BM[b], g_BM[b]);
    }
    SimTK_TEST_EQ_TOL(lhs, rhs, 1e-10);
}

int main() {
    SimTK_START_TEST("TestJacobians");
        SimTK_SUBTEST(testMultiplyByPositionJacobianWrtInboardFramePositions);
        SimTK_SUBTEST(
            testMultiplyByPositionJacobianWrtInboardFramePositionsTranspose);
        SimTK_SUBTEST(testMultiplyByPositionJacobianWrtOutboardFramePositions);
        SimTK_SUBTEST(
            testMultiplyByPositionJacobianWrtOutboardFramePositionsTranspose);
    SimTK_END_TEST();
}
