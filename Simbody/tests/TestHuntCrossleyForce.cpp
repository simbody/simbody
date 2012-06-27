/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
    const Real k1 = 1.0;
    const Real k2 = 2.0;
    const Real stiffness1 = std::pow(k1, 2.0/3.0);
    const Real stiffness2 = std::pow(k2, 2.0/3.0);
    const Real dissipation1 = 0.5;
    const Real dissipation2 = 1.0;
    const Real us1 = 1.0;
    const Real us2 = 0.7;
    const Real ud1 = 0.5;
    const Real ud2 = 0.2;
    const Real uv1 = 0.1;
    const Real uv2 = 0.05;
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Translation sphere(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, sphere, ContactGeometry::Sphere(radius), Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0))); // y < 0
    HuntCrossleyForce hc(forces, contacts, setIndex);
    hc.setBodyParameters(ContactSurfaceIndex(0), k1, dissipation1, us1, ud1, uv1);
    hc.setBodyParameters(ContactSurfaceIndex(1), k2, dissipation2, us2, ud2, uv2);
    const Real vt = 0.001;
    hc.setTransitionVelocity(vt);
    assertEqual(vt, hc.getTransitionVelocity());
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

    const Real us = 2.0*us1*us2/(us1+us2);
    const Real ud = 2.0*ud1*ud2/(ud1+ud2);
    const Real uv = 2.0*uv1*uv2/(uv1+uv2);
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = 0;
        if (depth > 0)
            fh = (4.0/3.0)*stiffness*depth*std::sqrt(radius*stiffness*depth);
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            const Real vrel = std::abs(v/vt);
            Real ff = (v < 0 ? 1 : -1)*fh*(std::min(vrel, 1.0)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*std::fabs(v));
            const Vec3 totalForce = gravity+Vec3(ff, fh, 0);
            expectedForce = SpatialVec(Vec3(0), Vec3(0));
            Vec3 contactPointInSphere = sphere.findStationAtGroundPoint(state, Vec3(0, -stiffness1*depth/(stiffness1+stiffness2), 0));
            sphere.applyForceToBodyPoint(state, contactPointInSphere, totalForce, expectedForce);
            SpatialVec actualForce = system.getRigidBodyForces(state, Stage::Dynamics)[sphere.getMobilizedBodyIndex()];
            assertEqual(actualForce[0], expectedForce[sphere.getMobilizedBodyIndex()][0]);
            assertEqual(actualForce[1], expectedForce[sphere.getMobilizedBodyIndex()][1]);
        }
    }
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
