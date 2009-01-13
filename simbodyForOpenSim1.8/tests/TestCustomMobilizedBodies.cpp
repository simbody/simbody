/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

template<>
void assertEqual(Vector val1, Vector val2) {
    ASSERT(val1.size() == val2.size());
    for (int i = 0; i < val1.size(); ++i)
        assertEqual(val1[i], val2[i]);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2) {
    assertEqual(val1[0], val2[0]);
    assertEqual(val1[1], val2[1]);
}

template<>
void assertEqual(Transform val1, Transform val2) {
    assertEqual(val1.T(), val2.T());
    ASSERT(val1.R().isSameRotationToWithinAngle(val2.R(), TOL));
}

/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Translation.
 */

class CustomTranslation : public MobilizedBody::Custom::Implementation {
public:
    CustomTranslation(SimbodyMatterSubsystem& matter) : Implementation(matter, 3, 3, 0) {
    }
    Implementation* clone() const {
        return new CustomTranslation(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const {
        ASSERT(nq == 3);
        return Transform(Vec3(q[0], q[1], q[2]));
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 3);
        return SpatialVec(Vec3(0), Vec3(u[0], u[1], u[2]));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 3);
        Vec3::updAs(f) = F[1];
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 3);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 3);
        Vec3::updAs(f) = Vec3(0);
    }
    void setQToFitTransform(const State&, const Transform& X_FM, int nq, Real* q) const {
        ASSERT(nq == 3);
        Vec3::updAs(q) = X_FM.T();
    }
    void setUToFitVelocity(const State&, const SpatialVec& V_FM, int nu, Real* u) const {
        ASSERT(nu == 3);
        Vec3::updAs(u) = V_FM[1];
    }
};

/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Ball.
 */

class CustomBall : public MobilizedBody::Custom::Implementation {
public:
    CustomBall(SimbodyMatterSubsystem& matter) : Implementation(matter, 3, 4, 4) {
    }
    Implementation* clone() const {
        return new CustomBall(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const {
        Transform t(Vec3(0));
        if (getUseEulerAngles(s)) {
            ASSERT(nq == 3);
            t.updR().setRotationToBodyFixedXYZ(Vec3::getAs(q)); 
        }
        else {
            ASSERT(nq == 4);
            t.updR().setRotationFromQuaternion(Quaternion(Vec4::getAs(q)));
        }
        return t;
        
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 3);
        return SpatialVec(Vec3(u[0], u[1], u[2]), Vec3(0));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 3);
        Vec3::updAs(f) = F[0];
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 3);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 3);
        Vec3::updAs(f) = Vec3(0);
    }
    void multiplyByQMatrix(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 3 && nOut == 3);
            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(&q[0])) * ~R_FM;
            if (transposeMatrix) Row3::updAs(out) = Row3::getAs(in) * Q;
            else                 Vec3::updAs(out) = Q * Vec3::getAs(in);
        }
        else {
            if (transposeMatrix)
                ASSERT(nIn == 4 && nOut == 3)
            else
                ASSERT(nIn == 3 && nOut == 4)
            const Mat43 Q = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) Row3::updAs(out) = Row4::getAs(in) * Q;
            else                 Vec4::updAs(out) = Q * Vec3::getAs(in);
        }
    }
    void multiplyByQInverse(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 3 && nOut == 3);
            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 QInv = R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(&q[0]));
            if (transposeMatrix) Row3::updAs(out) = Row3::getAs(in) * QInv;
            else                 Vec3::updAs(out) = QInv * Vec3::getAs(in);
        }
        else {
            if (transposeMatrix)
                ASSERT(nIn == 3 && nOut == 4)
            else
                ASSERT(nIn == 4 && nOut == 3)
            const Mat34 Q = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) Row4::updAs(out) = Row3::getAs(in) * Q;
            else                 Vec3::updAs(out) = Q * Vec4::getAs(in);
        }
    }
    void multiplyByQDotMatrix(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 3 && nOut == 3);
            const Rotation& R_FM = getMobilizerTransform(s).R();
            Vec3::updAs(out) = Rotation::convertAngVelDotToBodyFixed123DotDot(Vec3::getAs(&q[0]), ~R_FM*Vec3::getAs(in), Vec3(0));
        }
        else if (transposeMatrix) {
            ASSERT(nIn == 4 && nOut == 3)
            ASSERT(false); // I didn't bother to implement this case, since it currently is never used.
        }
        else {
            ASSERT(nIn == 3 && nOut == 4)
            Vec4::updAs(out) = Rotation::convertAngVelDotToQuaternionDotDot(Vec4::getAs(&q[0]), Vec3::getAs(in), Vec3(0));
        }
    }
    void setQToFitTransform(const State& s, const Transform& X_FM, int nq, Real* q) const {
        if (getUseEulerAngles(s)) {
            ASSERT(nq == 3);
            Vec3::updAs(q) = X_FM.R().convertRotationToBodyFixedXYZ();
        }
        else {
            ASSERT(nq == 4);
            Vec4::updAs(q) = X_FM.R().convertRotationToQuaternion().asVec4();
        }
    }
    void setUToFitVelocity(const State& s, const SpatialVec& V_FM, int nu, Real* u) const {
        ASSERT(nu == 3);
        Vec3::updAs(u) = V_FM[0];
    }
};

/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Free.
 */

class CustomFree : public MobilizedBody::Custom::Implementation {
public:
    CustomFree(SimbodyMatterSubsystem& matter) : Implementation(matter, 6, 7, 4) {
    }
    Implementation* clone() const {
        return new CustomFree(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const {
        Transform t(Vec3::getAs(&q[nq-3]));
        if (getUseEulerAngles(s)) {
            ASSERT(nq == 6);
            t.updR().setRotationToBodyFixedXYZ(Vec3::getAs(q)); 
        }
        else {
            ASSERT(nq == 7);
            t.updR().setRotationFromQuaternion(Quaternion(Vec4::getAs(q)));
        }
        return t;
        
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 6);
        return SpatialVec(Vec3(u[0], u[1], u[2]), Vec3(u[3], u[4], u[5]));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(f)) = F;
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const {
        ASSERT(nu == 6);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const {
        ASSERT(nu == 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(f)) = SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByQMatrix(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 6 && nOut == 6);
            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(&q[0])) * ~R_FM;
            if (transposeMatrix) Row3::updAs(out) = Row3::getAs(in) * Q;
            else                 Vec3::updAs(out) = Q * Vec3::getAs(in);
        }
        else {
            if (transposeMatrix)
                ASSERT(nIn == 7 && nOut == 6)
            else
                ASSERT(nIn == 6 && nOut == 7)
            const Mat43 Q = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) Row3::updAs(out) = Row4::getAs(in) * Q;
            else                 Vec4::updAs(out) = Q * Vec3::getAs(in);
        }
        Vec3::updAs(&out[nOut-3]) = Vec3::getAs(&in[nIn-3]);
   }
    void multiplyByQInverse(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 6 && nOut == 6);
            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 QInv = R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(&q[0]));
            if (transposeMatrix) Row3::updAs(out) = Row3::getAs(in) * QInv;
            else                 Vec3::updAs(out) = QInv * Vec3::getAs(in);
        }
        else {
            if (transposeMatrix)
                ASSERT(nIn == 6 && nOut == 7)
            else
                ASSERT(nIn == 7 && nOut == 6)
            const Mat34 Q = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) Row4::updAs(out) = Row3::getAs(in) * Q;
            else                 Vec3::updAs(out) = Q * Vec4::getAs(in);
        }
        Vec3::updAs(&out[nOut-3]) = Vec3::getAs(&in[nIn-3]);
    }
    void multiplyByQDotMatrix(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out) const {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            ASSERT(nIn == 6 && nOut == 6);
            const Rotation& R_FM = getMobilizerTransform(s).R();
            Vec3::updAs(out) = Rotation::convertAngVelDotToBodyFixed123DotDot(Vec3::getAs(&q[0]), ~R_FM*Vec3::getAs(in), Vec3(0));
        }
        else if (transposeMatrix) {
            ASSERT(nIn == 7 && nOut == 6)
            ASSERT(false); // I didn't bother to implement this case, since it currently is never used.
        }
        else {
            ASSERT(nIn == 6 && nOut == 7)
            Vec4::updAs(out) = Rotation::convertAngVelDotToQuaternionDotDot(Vec4::getAs(&q[0]), Vec3::getAs(in), Vec3(0));
        }
        Vec3::updAs(&out[nOut-3]) = Vec3(0);
    }
    void setQToFitTransform(const State& s, const Transform& X_FM, int nq, Real* q) const {
        if (getUseEulerAngles(s)) {
            ASSERT(nq == 6);
            Vec3::updAs(q) = X_FM.R().convertRotationToBodyFixedXYZ();
        }
        else {
            ASSERT(nq == 7);
            Vec4::updAs(q) = X_FM.R().convertRotationToQuaternion().asVec4();
        }
        Vec3::updAs(&q[nq-3]) = X_FM.T();
    }
    void setUToFitVelocity(const State& s, const SpatialVec& V_FM, int nu, Real* u) const {
        ASSERT(nu == 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(u)) = V_FM;
    }
};

void compareMobilizedBodies(const MobilizedBody& b1, const MobilizedBody& b2, bool eulerAngles, int expectedQ, int expectedU) {
    const SimbodyMatterSubsystem& matter = b1.getMatterSubsystem();
    const System& system = matter.getSystem();
    
    // Set whether to use Euler angles.
    
    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, eulerAngles);
    system.realizeModel(state);
    
    // Make sure the number of state variables is correct.
    
    assertEqual(b1.getNumQ(state), expectedQ);
    assertEqual(b1.getNumU(state), expectedU);
    assertEqual(b2.getNumQ(state), expectedQ);
    assertEqual(b2.getNumU(state), expectedU);

    // Set all the state variables to random values.

    Random::Gaussian random;
    int nq = state.getNQ()/2;
    for (int i = 0; i < nq; ++i)
        state.updQ()[i] = state.updQ()[i+nq] = random.getValue();
    int nu = state.getNU()/2;
    for (int i = 0; i < nu; ++i)
        state.updU()[i] = state.updU()[i+nu] = random.getValue();
    system.realize(state, Stage::Acceleration);
        
    // Compare state variables and their derivatives.
    
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        assertEqual(b1.getOneQ(state, i), b2.getOneQ(state, i));
        assertEqual(b1.getOneQDot(state, i), b2.getOneQDot(state, i));
        assertEqual(b1.getOneQDotDot(state, i), b2.getOneQDotDot(state, i));
    }
    for (int i = 0; i < b1.getNumU(state); ++i) {
        assertEqual(b1.getOneU(state, i), b2.getOneU(state, i));
        assertEqual(b1.getOneUDot(state, i), b2.getOneUDot(state, i));
    }
    
    // Compare lots of properties of the two bodies.
    
    assertEqual(b1.getBodyTransform(state), b2.getBodyTransform(state));
    assertEqual(b1.getBodyVelocity(state), b2.getBodyVelocity(state));
    assertEqual(b1.getBodyAcceleration(state), b2.getBodyAcceleration(state));
    assertEqual(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state));
    assertEqual(b1.getBodyOriginVelocity(state), b2.getBodyOriginVelocity(state));
    assertEqual(b1.getBodyOriginAcceleration(state), b2.getBodyOriginAcceleration(state));
    assertEqual(b1.getMobilizerTransform(state), b2.getMobilizerTransform(state));
    assertEqual(b1.getMobilizerVelocity(state), b2.getMobilizerVelocity(state));
    
    // Test methods that multiply by various matrices.
    
    Vector tempq(state.getNQ());
    Vector tempu(state.getNU());
    matter.multiplyByQMatrix(state, false, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    matter.multiplyByQMatrix(state, true, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByQMatrixInverse(state, false, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByQMatrixInverse(state, true, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    
    // Have them calculate q and u, and see if they agree.
    
    Transform t(Rotation(random.getValue(), Vec3(random.getValue(), random.getValue(), random.getValue())), Vec3(random.getValue(), random.getValue(), random.getValue()));
    b1.setQToFitTransform(state, t);
    b2.setQToFitTransform(state, t);
    assertEqual(b1.getQAsVector(state), b2.getQAsVector(state));
    SpatialVec v(Vec3(random.getValue(), random.getValue(), random.getValue()), Vec3(random.getValue(), random.getValue(), random.getValue()));
    b1.setUToFitVelocity(state, v);
    b2.setUToFitVelocity(state, v);
    assertEqual(b1.getUAsVector(state), b2.getUAsVector(state));
    
    // Simulate the system, and see if the two bodies remain identical.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    assertEqual(b1.getQAsVector(integ.getState()), b2.getQAsVector(integ.getState()));
    assertEqual(b1.getUAsVector(integ.getState()), b2.getUAsVector(integ.getState()));
}

void testCustomTranslation() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Translation t1(matter.Ground(), body);
    MobilizedBody::Translation t2(t1, body);
    MobilizedBody::Custom c1(matter.Ground(), new CustomTranslation(matter), body);
    MobilizedBody::Custom c2(c1, new CustomTranslation(matter), body);
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 3, 3);
    compareMobilizedBodies(t2, c2, true, 3, 3);
}

void testCustomBall() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Ball t1(matter.Ground(), body);
    MobilizedBody::Ball t2(t1, Vec3(1, 0, 0), body, Vec3(0, 2, 0));
    MobilizedBody::Custom c1(matter.Ground(), new CustomBall(matter), body);
    MobilizedBody::Custom c2(c1, new CustomBall(matter), Vec3(1, 0, 0), body, Vec3(0, 2, 0));
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 4, 3);
    compareMobilizedBodies(t2, c2, true, 3, 3);
}

void testCustomFree() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free t1(matter.Ground(), body);
    MobilizedBody::Free t2(t1, body);
    MobilizedBody::Custom c1(matter.Ground(), new CustomFree(matter), body);
    MobilizedBody::Custom c2(c1, new CustomFree(matter), body);
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 7, 6);
    compareMobilizedBodies(t2, c2, true, 6, 6);
}

int main() {
    try {
        testCustomTranslation();
        testCustomBall();
        testCustomFree();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
