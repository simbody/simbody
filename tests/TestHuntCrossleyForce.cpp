/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

void testForces() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    GeneralForceSubsystem forces(system);
    const Vec3 gravity = Vec3(0, -9.8, 0);
    Force::UniformGravity(forces, matter, gravity, 0);
    const Real radius = 0.8;
    const Real stiffness1 = 1.0;
    const Real stiffness2 = 2.0;
    const Real dissipation1 = 0.5;
    const Real dissipation2 = 1.0;
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Translation sphere(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, sphere, ContactGeometry::Sphere(radius), Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0))); // y < 0
    HuntCrossleyForce hc(forces, contacts, setIndex);
    hc.setBodyParameters(0, stiffness1, dissipation1);
    hc.setBodyParameters(1, stiffness2, dissipation2);
    State state = system.realizeTopology();
    
    // Position the sphere at a variety of positions and check the normal force.
    
    const Real stiffness = stiffness1*stiffness2/(stiffness1+stiffness2);
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        system.realize(state, Stage::Dynamics);
        const Real depth = radius-height;
        Real f = 0;
        if (depth > 0)
            f = (4.0/3.0)*stiffness*depth*std::sqrt(radius*stiffness*depth);
        assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()][1], gravity+Vec3(0, f, 0));
    }
    
    // Now do it with a vertical velocity and see if the dissipation force is correct.

    const Real dissipation = (dissipation1*stiffness2+dissipation2*stiffness1)/(stiffness1+stiffness2);
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = 0;
        if (depth > 0)
            fh = (4.0/3.0)*stiffness*depth*std::sqrt(radius*stiffness*depth);
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(0, -v, 0));
            system.realize(state, Stage::Dynamics);
            Real f = fh*(1.0+1.5*dissipation*v);
            if (f < 0)
                f = 0;
            assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()][1], gravity+Vec3(0, f, 0));
        }
    }
    
    // Do it with a horizontal velocity and see if the friction force is correct.
}

int main() {
    try {
        testForces();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
