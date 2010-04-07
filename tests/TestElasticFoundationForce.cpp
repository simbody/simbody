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

const Real TOL = 1e-5;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL || abs(val1-val2)/max(abs(val1), abs(val2)) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        assertEqual(val1[i], val2[i]);
}

void testForces() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    GeneralForceSubsystem forces(system);

    // Create a triangle mesh in the shape of a pyramid.
    
    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(1, 0, 1));
    vertices.push_back(Vec3(0, 0, 1));
    vertices.push_back(Vec3(0.5, 1, 0.5));
    vector<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {2, 1, 4}, {3, 2, 4}, {0, 3, 4}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);

    // Create the mobilized bodies and configure the contact model.
    
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, mesh, ContactGeometry::TriangleMesh(vertices, faceIndices), Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0))); // y < 0
    ElasticFoundationForce ef(forces, contacts, setIndex);
    Real stiffness = 1e9, dissipation = 0.01, us = 0.1, ud = 0.05, uv = 0.01, vt = 0.01;
    ef.setBodyParameters(0, stiffness, dissipation, us, ud, uv);
    ef.setTransitionVelocity(vt);
    ASSERT(ef.getTransitionVelocity() == vt);
    State state = system.realizeTopology();
    
    // Position the sphere at a variety of positions and check the normal force.
    
    for (Real depth = -0.1; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        system.realize(state, Stage::Dynamics);
        Real f = 0;
        if (depth > 0)
            f = stiffness*depth;
        assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()][1], Vec3(0, f, 0));
        assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[matter.getGround().getMobilizedBodyIndex()][1], Vec3(0, -f, 0));
    }
    
    // Now do it with a vertical velocity and see if the dissipation force is correct.

    for (Real depth = -0.105; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            mesh.setUToFitLinearVelocity(state, Vec3(0, -v, 0));
            system.realize(state, Stage::Dynamics);
            Real f = (depth > 0 ? stiffness*depth*(1+dissipation*v) : 0);
            if (f < 0)
                f = 0;
            assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()][1], Vec3(0, f, 0));
        }
    }
    
    // Do it with a horizontal velocity and see if the friction force is correct.

    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real depth = -0.105; depth < 0.1; depth += 0.01) {
        mesh.setQToFitTranslation(state, Vec3(0, -depth, 0));
        Real fh = 0;
        if (depth > 0)
            fh = stiffness*depth;
        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            mesh.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            const Real vrel = std::abs(v/vt);
            Real ff = (v < 0 ? 1 : -1)*fh*(std::min(vrel, 1.0)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*std::fabs(v));
            const Vec3 totalForce = Vec3(ff, fh, 0);
            expectedForce = SpatialVec(Vec3(0), Vec3(0));
            Vec3 contactPoint1 = mesh.findStationAtGroundPoint(state, Vec3(2.0/3.0, 0, 1.0/3.0));
            mesh.applyForceToBodyPoint(state, contactPoint1, 0.5*totalForce, expectedForce);
            Vec3 contactPoint2 = mesh.findStationAtGroundPoint(state, Vec3(1.0/3.0, 0, 2.0/3.0));
            mesh.applyForceToBodyPoint(state, contactPoint2, 0.5*totalForce, expectedForce);
            SpatialVec actualForce = system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()];
            assertEqual(actualForce[0], expectedForce[mesh.getMobilizedBodyIndex()][0]);
            assertEqual(actualForce[1], expectedForce[mesh.getMobilizedBodyIndex()][1]);
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
