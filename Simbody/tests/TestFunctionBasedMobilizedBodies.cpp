/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Ajay Seth                                          *
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

#include <vector>

using namespace SimTK;
using namespace std;

const Real TOL = 1e-6;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2, Real tol = TOL) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, Real tol = TOL) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template<>
void assertEqual(Vector val1, Vector val2, Real tol) {
    ASSERT(val1.size() == val2.size());
    for (int i = 0; i < val1.size(); ++i)
        assertEqual(val1[i], val2[i], tol);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2, Real tol) {
    assertEqual(val1[0], val2[0], tol);
    assertEqual(val1[1], val2[1], tol);
}

template<>
void assertEqual(Transform val1, Transform val2, Real tol) {
    assertEqual(val1.p(), val2.p(), tol);
    ASSERT(val1.R().isSameRotationToWithinAngle(val2.R(), tol));
}

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
        state.updU()[i] = state.updU()[i+nu] = (eulerAngles ? 0.0 : random.getValue());
    system.realize(state, Stage::Acceleration);
        
    // Compare state variables and their derivatives.
    
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        assertEqual(b1.getOneQ(state, i), b2.getOneQ(state, i));
        assertEqual(b1.getOneQDot(state, i), b2.getOneQDot(state, i));
        assertEqual(b1.getOneQDotDot(state, i), b2.getOneQDotDot(state, i));
    }
    /*
    for (int i = 0; i < b1.getNumU(state); ++i) {
        assertEqual(b1.getOneU(state, i), b2.getOneU(state, i));
        assertEqual(b1.getOneUDot(state, i), b2.getOneUDot(state, i));
    }
    */
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
    /*
    matter.multiplyByN(state, false, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    matter.multiplyByN(state, true, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByNInv(state, false, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByNInv(state, true, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    */
    
    // Have them calculate q and u, and see if they agree.
    
    if (!eulerAngles) { // The optimizer does not work reliably for Euler angles, since it can hit a singularity
        Transform t = b1.getBodyTransform(state);
        b1.setQFromVector(state, Vector(b1.getNumQ(state), 0.0));
        b2.setQFromVector(state, Vector(b2.getNumQ(state), 0.0));
        b1.setQToFitTransform(state, t);
        b2.setQToFitTransform(state, t);
        system.realize(state, Stage::Velocity);
        assertEqual(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state), 1e-2);
        assertEqual((~b1.getBodyRotation(state)*b2.getBodyRotation(state)).convertRotationToAngleAxis()[0], 0.0, 1e-2);
        SpatialVec v = b1.getBodyVelocity(state);
        b1.setUFromVector(state, Vector(b1.getNumU(state), 0.0));
        b2.setUFromVector(state, Vector(b2.getNumU(state), 0.0));
        b1.setUToFitVelocity(state, v);
        b2.setUToFitVelocity(state, v);
        assertEqual(b1.getUAsVector(state), b2.getUAsVector(state), 1e-2);
    }
    
    // Simulate the system, and see if the two bodies remain identical.
    b2.setQFromVector(state, b1.getQAsVector(state));
    b2.setUFromVector(state, b1.getUAsVector(state));
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    
    assertEqual(b1.getQAsVector(integ.getState()), b2.getQAsVector(integ.getState()));
    assertEqual(b1.getQDotAsVector(integ.getState()), b2.getQDotAsVector(integ.getState()));
}

class ConstantFunction : public Function {
// Implements a simple constant function, y = C
private:
    Real C;
public:

    //Default constructor
    ConstantFunction(){
        C = 0.0;
    }

    //Convenience constructor to specify constant value
    ConstantFunction(Real constant){
        C = constant;
    }

    Real calcValue(const Vector& x) const{
        return C;
    }

    // This is the pure virtual signature.
    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const{
        return 0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
    	return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const{
        // constant has no arguments
        return 0;
    }

    int getMaxDerivativeOrder() const{
        return 10;
    }
};

class LinearFunction : public Function {
// Implements a simple linear functional relationship, y = m*x + b
private:
    Real m;
    Real b;
public:

    //Default constructor
    LinearFunction(){
        m = 1.0;
        b = 0.0;
    }

    //Convenience constructor to specify the slope and Y-intercept of the linear r
    LinearFunction(Real slope, Real intercept){
        m = slope;
        b = intercept;
    }

    Real calcValue(const Vector& x) const{
        return m*x[0]+b;
    }

    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const{
        if (derivComponents.size() == 1)
            return m;
        return 0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
    	return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const{
        return 1;
    }

    int getMaxDerivativeOrder() const{
        return 10;
    }
};

class NonlinearFunction : public Function {
public:
    NonlinearFunction(){
    }
    Real calcValue(const Vector& x) const{
        return x[0]+x[1]*x[1];
    }
    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const{
        switch (derivComponents.size()) {
            case 1:
                return (derivComponents[0] == 0 ? 1.0 : x[1]);
            case 2:
                return (derivComponents[0] == 1 && derivComponents[1] == 1 ? 1.0 : 0.0);
        }
        return 0.0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
    	return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const{
        return 2;
    }
    int getMaxDerivativeOrder() const{
        return std::numeric_limits<int>::max();
    }
};

int defineMobilizerFunctions(std::vector<bool> &isdof, std::vector<std::vector<int> > &coordIndices, std::vector<const Function*> &functions1, std::vector<const Function*> &functions2)
{
    int nm = 0;
    for(int i=0; i<6; i++){
        if(isdof[i]) {
            std::vector<int> findex(1);
            findex[0] = nm++;
            functions1.push_back(new LinearFunction());
            functions2.push_back(new LinearFunction());
            coordIndices.push_back(findex);
        }
        else{
            std::vector<int> findex(0);
            functions1.push_back(new ConstantFunction());
            functions2.push_back(new ConstantFunction());
            coordIndices.push_back(findex);
        }
    }
    return nm;
}

void testFunctionBasedPin() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set the 1 spatial rotation about Z to be mobility
    isdof[2] = true;  //rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0.5), Inertia(1)));
    MobilizedBody::Pin p1(matter.Ground(), body);
    MobilizedBody::Pin p2(p1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

void testFunctionBasedSkewedPin() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);
    std::vector<Vec3> axes(6);

    // Set the 1 spatial rotation about first axis
    isdof[0] = true;  //rot 1
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    double angle = 0;

    axes[0] = Vec3(0,0,1);
    axes[1] = Vec3(0,1,0);
    axes[2] = Vec3(1,0,0);
    axes[3] = Vec3(1,0,0);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(0,0,1);

    Transform inParentPin = Transform(Rotation(angle, YAxis), Vec3(0));
    Transform inChildPin = Transform(Rotation(angle, YAxis), Vec3(0,1,0));
    
    Transform inParentFB = Transform(Vec3(0));
    Transform inChildFB = Transform(Vec3(0,1,0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    //Built-in
    MobilizedBody::Pin p1(matter.Ground(), inParentPin, body, inChildPin);
    MobilizedBody::Pin p2(p1, inParentPin, body, inChildPin);
    //Function-based
    MobilizedBody::FunctionBased fb1(matter.Ground(), inParentFB, body, inChildFB, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParentFB, body, inChildFB, nm, functions2, coordIndices, axes);

    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

void testFunctionBasedSlider() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set the 1 spatial translation along X to be mobility
    isdof[3] = true; //trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), body);
    MobilizedBody::Slider s2(s1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}


void testFunctionBasedSkewedSlider() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);
    std::vector<Vec3> axes(6);
    //axes[0] = Vec3(1/sqrt(2.0),0, -1/sqrt(2.0));
    axes[0] = Vec3(1,0,0);
    axes[1] = Vec3(0,1,0);
    axes[2] = Vec3(0,0,1);
    axes[3] = Vec3(0,0,1);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(1,0,0);
    // Set the 1 spatial translation along X to be mobility
    isdof[5] = true; //trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    Transform inParent = Transform(Vec3(0)); //Transform(Rotation(-Pi/2, YAxis));
    Transform inChild = Transform(Vec3(0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), inParent, body, inChild);
    MobilizedBody::Slider s2(s1, inParent, body, inChild);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices, axes);
    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}

void testFunctionBasedCylinder() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set 2 mobilities: rotation about and translation along Z
    isdof[2] = true;  //rot Z
    isdof[5] = true;  //trans Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Cylinder c1(matter.Ground(), body); 
    MobilizedBody::Cylinder c2(c1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(c2, fb2, false, nm, nm);
}

void testFunctionBasedUniversal() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set 2 rotation mobilities about body's X then Y
    isdof[0] = true;  //rot X
    isdof[1] = true;  //rot Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Universal u1(matter.Ground(), body); 
    MobilizedBody::Universal u2(u1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(u1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(u2, fb2, true, nm, nm);
}

void testFunctionBasedPlanar() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[2] = true;  //rot Z
    isdof[3] = true;  //trans X
    isdof[4] = true;  //trans Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    //Vec3(1/sqrt(2.0000000000000), 1/sqrt(2.0000000000000), 0)
    Body::Rigid body(MassProperties(1.0, Vec3(0, 0, 0), Inertia(1)));
    MobilizedBody::Planar u1(matter.Ground(), body); 
    MobilizedBody::Planar u2(u1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(u2, fb2, false, nm, nm);
}

void testFunctionBasedGimbal() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<bool> isdof(6,false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[0] = true;  //rot X
    isdof[1] = true;  //rot Y
    isdof[2] = true;  //rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::Gimbal b1(matter.Ground(), body); 
    MobilizedBody::Gimbal b2(b1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(b2, fb2, true, nm, nm);
}

void testFunctionBasedGimbalUserAxes() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;
    std::vector<Vec3> axes(6);
    std::vector<bool> isdof(6,false);

    isdof[0] = true;  //rot 1
    isdof[1] = true;  //rot 2
    isdof[2] = true;  //rot 3
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    // Sherm 20130213: I replaced the random number generator with some
    // firm numbers to prevent singularities from occurring on some platforms
    // based on different random number output.

    axes[0] = Vec3(0.05, 1.4,  0);
    axes[1] = Vec3(0.6,   0, -1.2);
    axes[2] = Vec3(0,     2, -0.055);
    axes[3] = Vec3(1,0,0);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(0,0,1);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    //Use massless bodies for generationg skewed-axes
    Body::Massless massLessBody;

    Transform inParent = Transform(Vec3(0));
    Transform inChild = Transform(Vec3(0,1,0));

       // Compared to standard built-in pin mobilizers with skewed axes
    // Pin rotates about Z-axis and need to align with first axis
    Transform parentPinAxis0 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis0 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0,0,0));
    Transform parentPinAxis1 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis1 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0,0,0));
    Transform parentPinAxis2 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis2 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0,1,0));
    
    MobilizedBody::Pin masslessPin0(matter.Ground(), parentPinAxis0, massLessBody, childPinAxis0);
    MobilizedBody::Pin masslessPin1(masslessPin0, parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin b1(masslessPin1, parentPinAxis2, body, childPinAxis2);

    MobilizedBody::Pin masslessPin00(b1, parentPinAxis0, massLessBody, childPinAxis0);
    MobilizedBody::Pin masslessPin01(masslessPin00, parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin b2(masslessPin01, parentPinAxis2, body, childPinAxis2);

    MobilizedBody::FunctionBased fb1(matter.Ground(), inParent, body, inChild, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParent, body, inChild, nm, functions2, coordIndices, axes);
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    // These were generated randomly but we want repeatability across 
    // machines so we'll use the same numbers every time. Note that we'll
    // re-use each of these twice, once for the pin joint system and once
    // for the function based mobilizers.
    Real initq[] = {1.41292, 0.048025,-1.19474, 0.618909,-0.0552235, 2.043930};
    Real initu[] = {1.53485, 0.546119,-1.55779,-1.872230, 0.0982929, 0.118798};
    int nq = state.getNQ()/2;
    assert(nq <= 6); // make sure we have enough random numbers!
    for (int i = 0; i < nq; ++i)
        state.updQ()[i] = state.updQ()[i+nq] = initq[i];
    int nu = state.getNU()/2;
    for (int i = 0; i < nu; ++i)
        state.updU()[i] = state.updU()[i+nu] = initu[i]; 

    system.realize(state, Stage::Acceleration);

    Transform Xb2 = b2.getBodyTransform(state);
    Transform Xfb2 = fb2.getBodyTransform(state);

    assertEqual(Xb2, Xfb2);

    SpatialVec A_b2 = b2.getBodyAcceleration(state);
    SpatialVec A_fb2 = fb2.getBodyAcceleration(state);

    assertEqual(A_b2, A_fb2);
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    const State& result = ts.getState();

    Vec3 com_bin = b2.getBodyOriginLocation(result);
    Vec3 com_fb = fb2.getBodyOriginLocation(result);
    
    assertEqual(com_bin, com_fb);
    assertEqual(b2.getBodyVelocity(result), fb2.getBodyVelocity(result));
    
    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    assertEqual(b2.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}

void testFunctionBasedTranslation() {
    // Test against built-in Translation mobilizer
    // for a total of 3 coordinates and 3 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof(6,true);
    //No rotations
    isdof[0] = false;  //rot X
    isdof[1] = false;  //rot Y
    isdof[2] = false;  //rot Z

    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    //Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));

    //Built-in mobilized bodies
    MobilizedBody::Translation b1(matter.Ground(), body);  
    MobilizedBody::Translation b2(b1, body);
    
    // Function-based
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    int nq = state.getNQ()/2;
    for (int i = 0; i < nq; ++i)
        state.updQ()[i] = state.updQ()[i+nq] = random.getValue();
    int nu = state.getNU()/2;
    for (int i = 0; i < nu; ++i)
        state.updU()[i] = state.updU()[i+nu] = random.getValue(); //0.0; //

    system.realize(state, Stage::Acceleration);

    Transform Xb2 = b2.getBodyTransform(state);
    Transform Xfb2 = fb2.getBodyTransform(state);

    SpatialVec A_b2 = b2.getBodyAcceleration(state);
    SpatialVec A_fb2 = fb2.getBodyAcceleration(state);

    assertEqual(A_b2, A_fb2);
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    const State& result = ts.getState();

    Vec3 com_bin = b2.getBodyOriginLocation(result);
    Vec3 com_fb = fb2.getBodyOriginLocation(result);
    
    assertEqual(com_bin, com_fb);
    assertEqual(b2.getBodyVelocity(result), fb2.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    assertEqual(b2.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}


void testFunctionBasedFree() {
    // Test against free joint using Euler angles for orientation (q)
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function*> functions1, functions2;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof(6,true);

    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0.2, -0.5, 0.1), Inertia(1.2)));

    //Built-in free
    MobilizedBody::Free b1(matter.Ground(), body);

    //Function-based equivalent?    
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
 
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    int nq = state.getNQ();
    nq = nq-nm;
  
    assert(nm == state.getNU()/2);

    // Get random q's and u's and set equivalent on both bodies.
    // (Not really random so we can get repeatability on all platforms.)
    Real initq[] = {0.455189,-0.383271,1.21353,-0.510623,-1.71438,0.968387};
    assert(nm <= 6);
    for (int i = 0; i < nm; ++i){
        // Free has slots for 4 rot q's and fb only has 3
        state.updQ()[i] = state.updQ()[i+nq] = initq[i]; //
    }

    system.realize(state, Stage::Position);
    SpatialVec inputVelocity(Vec3(-0.962157,0.523767,1.94993),
                             Vec3(-1.15752,0.436991,-0.787116));

    b1.setUToFitVelocity(state, inputVelocity);
    fb1.setUToFitVelocity(state, inputVelocity);

    system.realize(state, Stage::Acceleration);

    //cout << system.getRigidBodyForces(state, Stage::Dynamics)[b1.getMobilizedBodyIndex()] << endl;
    //cout << system.getRigidBodyForces(state, Stage::Dynamics)[fb1.getMobilizedBodyIndex()] << endl;

    Transform Xb1 = b1.getBodyTransform(state);
    Transform Xfb1 = fb1.getBodyTransform(state);

    assertEqual(Xb1, Xfb1); 

    SpatialVec A_b1 = b1.getBodyAcceleration(state);
    SpatialVec A_fb1 = fb1.getBodyAcceleration(state);

    assertEqual(A_b1, A_fb1);
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    const State &result = ts.getState();

    //cout << "Free and Function-based Us" << endl;
    //cout << result.getU()(0, 6) << endl;
    //cout << result.getU()(6, 6) << endl;

    Xb1 = b1.getBodyTransform(result);
    Xfb1 = fb1.getBodyTransform(result);

    Vec3 com_bin = b1.getBodyOriginLocation(result);
    Vec3 com_fb = fb1.getBodyOriginLocation(result);
    
    assertEqual(com_bin, com_fb);
    assertEqual(b1.getBodyVelocity(result), fb1.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    assertEqual(b1.getBodyAcceleration(result), fb1.getBodyAcceleration(result));
}

void testFunctionBasedFreeVsTranslationGimbal() {
    // Test function-based free against a combination of Translation and Gimbal mobilizer
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices1, coordIndices2a, coordIndices2b;
    std::vector<const Function*> functions1, temp;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof1(6,true), isdof2a(6,true), isdof2b(6,true);


    int nm1 = defineMobilizerFunctions(isdof1, coordIndices1, functions1, temp);
    int nm2a = 3;
    int nm2b = 3;

    // Check that we constructed the correct number of functions
    assert(nm1 == nm2a+nm2b);

    //Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(0.1, Vec3(0.25, -0.5, 0.1), Inertia(0.5)));
   
    // One Free-like function-based mmobilizer
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm1, functions1, coordIndices1);

    // Two function-based mmobilizers: 2a for translation and 2b for rotation 
    MobilizedBody::Translation massLessTrans(matter.Ground(), massLessBody);
    MobilizedBody::Gimbal b1(massLessTrans, body);

    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    int nq = state.getNQ();
    nq = nq/2;
    assert(nq == nm1);

    // Set rotation states first
    for (int i = 0; i < nm2b; ++i){
        state.updQ()[i] = state.updQ()[i+nq+nm2a] = 0.0; //random.getValue();
        state.updU()[i] = state.updU()[i+nq+nm2a] = random.getValue(); 
    }

    // Set translations states second
    for (int i = 0; i < nm2a; ++i){
        state.updQ()[i+nm2a] = state.updQ()[i+nq] = random.getValue();
        state.updU()[i+nm2a] = state.updU()[i+nq] = random.getValue(); 
    }

    system.realize(state, Stage::Acceleration);

    Transform Xfb1 = fb1.getBodyTransform(state);
    Transform Xb1 = b1.getBodyTransform(state);

    assertEqual(Xfb1, Xb1);

    SpatialVec A_fb1 = fb1.getBodyAcceleration(state);
    SpatialVec A_b1 = b1.getBodyAcceleration(state);

    assertEqual(A_fb1, A_b1);
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);

    const State& result = ts.getState();

    Xfb1 = fb1.getBodyTransform(result);
    Xb1 = b1.getBodyTransform(result);
    
    assertEqual(Xfb1, Xb1);
    assertEqual(fb1.getBodyVelocity(result), b1.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    assertEqual(fb1.getBodyAcceleration(result), b1.getBodyAcceleration(result));
}


void testFunctionBasedFreeVs2FunctionBased() {
    // Test against free joint that is a combination of Translation and Gimbal mobilizer
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices1, coordIndices2a, coordIndices2b;
    std::vector<const Function*> functions1, temp;
    std::vector<const Function*> functions2a, functions2b;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof1(6,true), isdof2a(6,true), isdof2b(6,true);

    // Just translation 
    isdof2a[0] = false;  //rot X
    isdof2a[1] = false;  //rot Y
    isdof2a[2] = false;  //rot Z

    // Just rotation
    isdof2b[3] = false;  //trans X
    isdof2b[4] = false;  //trans Y
    isdof2b[5] = false;  //trans Z

    int nm1 = defineMobilizerFunctions(isdof1, coordIndices1, functions1, temp);
    int nm2a = defineMobilizerFunctions(isdof2a, coordIndices2a, functions2a, temp);
    int nm2b = defineMobilizerFunctions(isdof2b, coordIndices2b, functions2b, temp);

    // Check that we constructed the correct number of functions
    assert(nm1 == nm2a+nm2b);

    //Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(0.1, Vec3(0.25, -0.5, 0.1), Inertia(0.5)));
   
    // One Free-like function-based mmobilizer
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm1, functions1, coordIndices1);

    // Two function-based mmobilizers: 2a for translation and 2b for rotation 
    MobilizedBody::FunctionBased massLessfb(matter.Ground(), massLessBody, nm2a, functions2a, coordIndices2a);
    MobilizedBody::FunctionBased fb2(massLessfb, body, nm2b, functions2b, coordIndices2b);

    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    int nq = state.getNQ();
    nq = nq/2;
    assert(nq == nm1);

    // Set rotation states first
    for (int i = 0; i < nm2b; ++i){
        state.updQ()[i] = state.updQ()[i+nq+nm2a] = random.getValue();
        state.updU()[i] = state.updU()[i+nq+nm2a] = random.getValue(); 
    }

    // Set translations states second
    for (int i = 0; i < nm2a; ++i){
        state.updQ()[i+nm2a] = state.updQ()[i+nq] = random.getValue();
        state.updU()[i+nm2a] = state.updU()[i+nq] = random.getValue(); 
    }

    system.realize(state, Stage::Acceleration);

    Transform Xfb1 = fb1.getBodyTransform(state);
    Transform Xfb2 = fb2.getBodyTransform(state);

    SpatialVec A_fb1 = fb1.getBodyAcceleration(state);
    SpatialVec A_fb2 = fb2.getBodyAcceleration(state);

    assertEqual(A_fb1, A_fb2);
    
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);

    const State& result = ts.getState();

    Vec3 com_fb1 = fb1.getBodyOriginLocation(result);
    Vec3 com_fb2 = fb2.getBodyOriginLocation(result);
    
    assertEqual(com_fb1, com_fb2);
    assertEqual(fb1.getBodyVelocity(result), fb2.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    assertEqual(fb1.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}



/**
 * Test a mobilized body based on functions that take multiple arguments.
 */
void testMultipleArguments() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    
    std::vector<std::vector<int> > coordIndices(6);
    std::vector<const Function*> functions(6);
    Vector coeff(3);
    coeff[0] = 0.5;
    coeff[1] = -0.5;
    coeff[2] = 1.0;
    functions[0] = new Function::Constant(0.0, 0);
    functions[1] = new Function::Constant(0.0, 0);
    functions[2] = new Function::Constant(0.0, 0);
    functions[3] = new NonlinearFunction();
    functions[4] = new Function::Linear(coeff);
    functions[5] = new Function::Constant(0.0, 0);
    coordIndices[3].push_back(0);
    coordIndices[3].push_back(1);
    coordIndices[4].push_back(0);
    coordIndices[4].push_back(1);

    // Create the system.

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::FunctionBased fb(matter.Ground(), body, 2, functions, coordIndices);
    State state = system.realizeTopology();
    
    // See if coordinates and velocities are calculated correctly.
    
    ASSERT(state.getNQ() == 2);
    ASSERT(state.getNU() == 2);
    state.updQ()[0] = 2.0;
    state.updQ()[1] = -3.0;
    state.updU()[0] = 0.1;
    state.updU()[1] = -0.4;
    system.realize(state, Stage::Acceleration);
    assertEqual(fb.getBodyTransform(state), Transform(Vec3(11.0, 3.5, 0.0)));
    assertEqual(fb.getBodyVelocity(state), SpatialVec(Vec3(0.0), Vec3(0.1+3.0*0.4, 0.25, 0.0)));

    // Simulate it.
    
    Real energy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5.0);
    assertEqual(energy, system.calcEnergy(ts.getState()));
}

int main() {
    try {
        cout << "FunctionBased MobilizedBodies vs. Built-in Types: " << endl;
      
        testFunctionBasedPin();
        cout << "Pin: Passed" << endl;

        testFunctionBasedSkewedPin();
        cout << "Skewed Pin: Passed" << endl;

        testFunctionBasedSlider();
        cout << "Slider: Passed" << endl;

        testFunctionBasedSkewedSlider();
        cout << "Skewed Slider: Passed" << endl;
        
        testFunctionBasedCylinder();
        cout << "Cylinder: Passed" << endl;

        testFunctionBasedUniversal();
        cout << "Universal: Passed" << endl;

        testFunctionBasedPlanar();
        cout << "Planar: Passed" << endl;

        testFunctionBasedGimbal();
        cout << "Gimbal: Passed" << endl;

        testMultipleArguments();
        cout << "Functions with Multiple Arguments: Passed" << endl;

        testFunctionBasedGimbalUserAxes();
        cout << "Gimbal with User Axes: Passed" << endl;

        testFunctionBasedTranslation();
        cout << "Translation: Passed" << endl;

        testFunctionBasedFree();
        cout << "Free: Passed" << endl;

        testFunctionBasedFreeVsTranslationGimbal();
        cout << "F-B Free vs. Translation and Gimbal: Passed" << endl;

        testFunctionBasedFreeVs2FunctionBased();
        cout << "F-B Free vs. Combination of Two Function-Based: Passed" << endl;
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
