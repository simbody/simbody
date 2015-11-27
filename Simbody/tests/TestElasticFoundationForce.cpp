/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Guillaume Jacquenot                                *
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

    // Create a triangle mesh in the shape of a pyramid, with the
    // square base having area 1 (split into two triangles).

    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(1, 0, 1));
    vertices.push_back(Vec3(0, 0, 1));
    vertices.push_back(Vec3(0.5, 1, 0.5));
    vector<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4},
                       {2, 1, 4}, {3, 2, 4}, {0, 3, 4}};
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
    ef.setBodyParameters(ContactSurfaceIndex(0), stiffness, dissipation, us, ud, uv);
    ef.setTransitionVelocity(vt);
    ASSERT(ef.getTransitionVelocity() == vt);
    State state = system.realizeTopology();

    // Position the pyramid at a variety of positions and check the normal
    // force.

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

/**
 * @brief This test compares the numerical result of a sphere
 *        in contact with a plane, using the elastic foundation
 *        model.
 *        The analytical solution of this problem is given by
 *        the product of the stiffness with the volume of
 *        the sphere in the plane, i.e the volume of a spherical
 *        cap.
 *        The volume of a spherical cap is:
 *          Vcap = Pi*h*h/3.0*(3.0*r-h)
 *        where
 *          r is the radis of the sphere
 *          h the height of the cap. In our case, the penetration
 *          depth
 * @note If we want to go further, we can observe that doubling
 *       the penetration depth results in multiplying the normal
 *       effort by 4.
 *       This is different from Hertz theory, where doubling the
 *       penetration depth results in multiplying
 *       the normal effort by 2^(3/2)~2.68
 *
 */
void testEffSphereOnPlaneOldFormulation(bool verbose = false)
{
    // Material properties for sphere
    const Real stiffness = 1e9;
    const Real dissipation = 0.0, us = 0.0, ud = 0.0, uv = 0.0, vt = 0.0;
    // Sphere radius
    const Real radius = 1.0;
    // Define initial penetration
    const Real initialPenetration = 0.002;
    // Define the number of tests to perform
    const int maxLevel = 6;
    // Define some tolerances for each level in %
    const Real tolerances[6]= {0.15, 0.07, 0.03, 0.02, 0.01, 0.02};
    for (int i=0;i<maxLevel;++i)
    {
        // For each level, penetration is double
        const Real penetration = initialPenetration * pow(2.0,(Real)i);
        // Creation of the classical problem
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        GeneralContactSubsystem contacts(system);
        GeneralForceSubsystem forces(system);
        const ContactSetIndex setIndex = contacts.createContactSet();
        // Creation a sphere with 6 levels of refinement
        const PolygonalMesh sphereMesh(PolygonalMesh::createSphereMesh(radius, 6));
        // Create the mobilized bodies and configure the contact model.
        const Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
        const MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());
        contacts.addBody(setIndex, mesh, ContactGeometry::TriangleMesh(sphereMesh), Transform());
        contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(),
                         Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0.0, penetration-radius, 0.0))); // y < penetration-radius
        ElasticFoundationForce ef(forces, contacts, setIndex);
        ef.setBodyParameters(ContactSurfaceIndex(0), stiffness, dissipation, us, ud, uv);
        ef.setTransitionVelocity(vt);
        const State state = system.realizeTopology();
        system.realize(state, Stage::Dynamics);
        const SpatialVec r = system.getRigidBodyForces(state, Stage::Dynamics)[mesh.getMobilizedBodyIndex()];
        const Real volumeSphericalCap = Pi*penetration*penetration/3.0*(3.0*radius-penetration);
        const Real theoreticalResult = stiffness*volumeSphericalCap;
        const Real numericalResult = r[1][1];
        ASSERT(abs(r[1][0])<TOL);
        ASSERT(abs(r[1][2])<TOL);
        const Real relativeDifference = abs((numericalResult/theoreticalResult)-1.0);
        if (verbose) {
            cout<<"Effort for penetration : "
                <<penetration*1000.0<<" mm -> F = "<<numericalResult<<" N "
                <<"(theoretical result : "<<theoreticalResult<< " N "
                <<" relative difference : "<<100.0*relativeDifference<<" %)"<<endl;
        }
        ASSERT(abs((numericalResult/theoreticalResult)-1.0)<tolerances[i]);
    }
}

void testEffSphereOnPlaneNewFormulation(bool verbose = false)
{
    // Global stiffness of the contact: each material will have
    // twice this stiffness to obtain this global stiffness in the contact
    // 1/kG = 1/k1 + 1/k2
    const Real stiffness = 1e9;
    const Real dissipation = 0.0, us = 0.0, ud = 0.0, uv = 0.0;
    const Real vt = 1.0e-2;
    // Sphere radius
    const Real radius = 1.0;
    // Define initial penetration
    const Real initialPenetration = 0.002;
    const int maxLevel = 6;
    // Define some tolerances for each level in %
    const Real tolerances[6]= {0.15, 0.07, 0.03, 0.02, 0.01, 0.02};
    for (int i=0;i<maxLevel;++i)
    {
        // For each level, penetration is double
        const Real penetration = initialPenetration * pow(2.0,(Real)i);
        // Creation of the classical problem
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        ContactTrackerSubsystem tracker(system);
        CompliantContactSubsystem contactForces(system, tracker);
        contactForces.setTransitionVelocity(vt);
        matter.Ground().updBody().addContactSurface(
            Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0.0,penetration-radius,0.0)), // y < penetration-radius
            ContactSurface(ContactGeometry::HalfSpace(),
                           ContactMaterial(2.0*stiffness, dissipation, us, ud, uv),
                           1.0));
        Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
        body.addContactSurface(Transform(),
            ContactSurface(ContactGeometry::TriangleMesh(PolygonalMesh::createSphereMesh(radius, 6)),
                           ContactMaterial(2.0*stiffness, dissipation, us, ud, uv),
                           1.0));
        const MobilizedBody::Translation mesh(matter.updGround(), Transform(), body, Transform());
        const State state = system.realizeTopology();
        system.realize(state, Stage::Dynamics);
        if (verbose) {
            cout << "Num contacts: " << contactForces.getNumContactForces(state) << endl;
        }
        ASSERT(contactForces.getNumContactForces(state)==1);
        const ContactForce& force = contactForces.getContactForce(state,0);
        const Vec3& frc = force.getForceOnSurface2()[1];
        ASSERT(abs(frc[0])<TOL);
        ASSERT(abs(frc[2])<TOL);
        const Real numericalResult = frc[1];
        const Real volumeSphericalCap = Pi*penetration*penetration/3.0*(3.0*radius-penetration);
        const Real theoreticalResult = stiffness*volumeSphericalCap;
        const Real relativeDifference = abs((numericalResult/theoreticalResult)-1.0);
        if (verbose) {
            cout<<force;
            cout<<"Effort for penetration : "
                <<penetration*1000.0<<" mm -> F = "<<numericalResult<<" N "
                <<"(theoretical result : "<<theoreticalResult<< " N "
                <<" relative difference : "<<100.0*relativeDifference<<" %)"<<endl;
        }
        ASSERT(abs((numericalResult/theoreticalResult)-1.0)<tolerances[i]);
    }
}

int main() {
    try {
        testForces();
        testEffSphereOnPlaneOldFormulation();
        testEffSphereOnPlaneNewFormulation();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
