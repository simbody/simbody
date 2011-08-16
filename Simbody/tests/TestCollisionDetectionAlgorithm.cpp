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

#include <set>

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

void verifyPointContact(const Array_<Contact>& contacts, int surface1, int surface2, const Vec3& normal, const Vec3& location, Real depth, Real r1, Real r2) {
    ASSERT(contacts.size() == 1);
    ASSERT(PointContact::isInstance(contacts[0]));
    const PointContact& c = static_cast<const PointContact&>(contacts[0]);
    assertEqual((int) c.getSurface1(), surface1);
    assertEqual((int) c.getSurface2(), surface2);
    assertEqual(c.getNormal(), normal);
    assertEqual(c.getDepth(), depth);
    assertEqual(min(c.getRadiusOfCurvature1(), c.getRadiusOfCurvature2()), min(r1, r2));
    assertEqual(max(c.getRadiusOfCurvature1(), c.getRadiusOfCurvature2()), max(r1, r2));
    assertEqual(c.getEffectiveRadiusOfCurvature(), sqrt(r1*r2));
    assertEqual(c.getLocation(), location);
}

void testHalfSpaceSphere() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Real radius = 0.8;
    Vec3 center(0.1, -0.3, 0.3);
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free sphere(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, sphere, ContactGeometry::Sphere(radius), center);
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();
    Vec3 centerInGround;
    for (int iteration = 0; iteration < 100; ++iteration) {
        // Pick a random positions for the sphere.

        for (int i = 0; i < state.getNY(); i++)
            state.updY()[i] = 5*random.getValue();
        system.realize(state, Stage::Dynamics);
        centerInGround = sphere.findStationLocationInGround(state, center);
        
        // Check the results of collision detection.
        
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);
        if (centerInGround[1] > radius+1) {
            ASSERT(contact.size() == 0);
        }
        else {
            Real depth = radius-centerInGround[1]+1;
            verifyPointContact(contact, 1, 0, Vec3(0, 1, 0), Vec3(centerInGround[0], 1-0.5*depth, centerInGround[2]), depth, radius, radius);
        }
    }
}

void testSphereSphere() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    const int numBodies = 10;
    Real radius[numBodies];
    Vec3 center[numBodies];
    Random::Uniform random(0.0, 1.0);
    for (int i = 0; i < numBodies; i++) {
        radius[i] = random.getValue();
        center[i] = Vec3(random.getValue(), random.getValue(), random.getValue());
    }
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    for (int i = 0; i < numBodies; ++i) {
        MobilizedBody::Free b(matter.updGround(), Transform(), body, Transform());
        contacts.addBody(setIndex, b, ContactGeometry::Sphere(radius[i]), center[i]);
    }
    State state = system.realizeTopology();
    Vec3 centerInGround[numBodies];
    for (int iteration = 0; iteration < 100; ++iteration) {
        // Pick random positions for all the bodies.

        for (int i = 0; i < state.getNY(); i++)
            state.updY()[i] = 5*random.getValue();
        system.realize(state, Stage::Dynamics);
        for (MobilizedBodyIndex index(1); index <= numBodies; ++index)
            centerInGround[index-1] = matter.getMobilizedBody(index).findStationLocationInGround(state, center[index-1]);
        
        // Make sure all contacts are accurate.
        
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);
        for (int i = 0; i < (int) contact.size(); i++) {
            ASSERT(PointContact::isInstance(contact[i]));
            const PointContact& c = static_cast<const PointContact&>(contact[i]);
            int body1 = c.getSurface1();
            int body2 = c.getSurface2();
            Vec3 delta = centerInGround[body2]-centerInGround[body1];
            assertEqual(delta.normalize(), c.getNormal());
            assertEqual(delta.norm(), radius[body1]+radius[body2]-c.getDepth());
            double r = radius[body1]*radius[body2]/(radius[body1]+radius[body2]);
            assertEqual(r, c.getRadiusOfCurvature1());
            assertEqual(r, c.getRadiusOfCurvature2());
        }

        // Make sure no contacts were missed.
        
        int expectedContacts = 0;
        for (int i = 0; i < numBodies; i++)
            for (int j = 0; j < i; j++)
                if ((centerInGround[i]-centerInGround[j]).norm() < radius[i]+radius[j])
                    expectedContacts++;
        ASSERT(contact.size() == expectedContacts);
    }
}

void testHalfSpaceEllipsoid() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Vec3 radii(0.8, 1.5, 2.1);
    Vec3 center(0.1, -0.3, 0.3); // Major axes span the ranges [-0.7, 0.9], [-1.8, 1.2], [-1.8, 2.4]
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free ellipsoid(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, ellipsoid, ContactGeometry::Ellipsoid(radii), center);
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();

    // Test a variety of positions.

    ellipsoid.setQToFitTransform(state, Transform(Rotation(), Vec3(0, 2.9, 0))); // [-0.7, 0.9], [1.1, 4.1], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    ellipsoid.setQToFitTransform(state, Transform(Rotation(), Vec3(0, 2.6, 0))); // [-0.7, 0.9], [0.8, 3.8], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    verifyPointContact(contacts.getContacts(state, setIndex), 1, 0, Vec3(0, 1, 0), Vec3(0.1, 0.9, 0.3), 0.2, 0.8, 2.1);
    ellipsoid.setQToFitTransform(state, Transform(Rotation(SimTK_PI/2, ZAxis), Vec3(0, 1.6, 0))); // [-1.2, 1.8], [0.9, 2.5], [-1.8, 2.4]
    system.realize(state, Stage::Dynamics);
    verifyPointContact(contacts.getContacts(state, setIndex), 1, 0, Vec3(0, 1, 0), Vec3(0.3, 0.95, 0.3), 0.1, 1.5, 2.1);
    ellipsoid.setQToFitTransform(state, Transform(Rotation(SimTK_PI/4, XAxis), Vec3(0, 3.1, 0)));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 1);
    const PointContact& c = static_cast<const PointContact&>(contacts.getContacts(state, setIndex)[0]);
    assertEqual(c.getNormal(), Vec3(0, 1, 0));
    ASSERT(c.getDepth() < 0.3);
    assertEqual(min(c.getRadiusOfCurvature1(), c.getRadiusOfCurvature2()), 0.8);
    ASSERT(max(c.getRadiusOfCurvature1(), c.getRadiusOfCurvature2()) > 1.5 && max(c.getRadiusOfCurvature1(), c.getRadiusOfCurvature2()) < 2.1);
    assertEqual(c.getLocation()[0], 0.1);
    assertEqual(c.getLocation()[1], 1-c.getDepth()/2);
    ASSERT(c.getLocation()[2] > 0);
}

bool verifyEllipsoidContact(const Contact& contact, const Vec3& radii1, const Vec3& radii2, const Vec3& center1, const Vec3& center2, const Transform& t1, const Transform& t2) {
    ASSERT(PointContact::isInstance(contact));
    const PointContact& c = static_cast<const PointContact&>(contact);

    // The "contact point" should be midway between the two surfaces along the normal direction.  Verify that.

    Vec3 loc1 = ~t1*(c.getLocation()+0.5*c.getDepth()*c.getNormal())-center1;
    assertEqual(loc1[0]*loc1[0]/(radii1[0]*radii1[0])+loc1[1]*loc1[1]/(radii1[1]*radii1[1])+loc1[2]*loc1[2]/(radii1[2]*radii1[2]), 1.0);
    Vec3 loc2 = ~t2*(c.getLocation()-0.5*c.getDepth()*c.getNormal())-center2;
    assertEqual(loc2[0]*loc2[0]/(radii2[0]*radii2[0])+loc2[1]*loc2[1]/(radii2[1]*radii2[1])+loc2[2]*loc2[2]/(radii2[2]*radii2[2]), 1.0);

    // Check that the normals are correct.  This test may occassionally fail (when points of very high
    // curvate cause the Newton iteration not to converge), so instead of an assertion, which just return
    // whether the normals were correct.

    UnitVec3 norm1(loc1[0]/(radii1[0]*radii1[0]), loc1[1]/(radii1[1]*radii1[1]), loc1[2]/(radii1[2]*radii1[2]));
    if (~norm1*(~t1.R()*c.getNormal()) < 0.999)
        return false;
    UnitVec3 norm2(loc2[0]/(radii2[0]*radii2[0]), loc2[1]/(radii2[1]*radii2[1]), loc2[2]/(radii2[2]*radii2[2]));
    if (-~norm2*(~t2.R()*c.getNormal()) < 0.999)
        return false;
    return true;
}

void testEllipsoidEllipsoid() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Vec3 radii1(0.8, 1.5, 2.1);
    Vec3 radii2(1.0, 1.2, 1.4);
    Vec3 center1(0, -0.2, 0.5); // Major axes span the ranges [-0.8, 0.8], [-1.7, 1.3], [-1.6, 2.6]
    Vec3 center2(0.1, 0, 0.3); // Major axes span the ranges [-0.9, 1.1], [-1.2, 1.2], [-1.1, 1.7]
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    Body::Rigid body2(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free ellipsoid1(matter.updGround(), Transform(), body, Transform());
    MobilizedBody::Free ellipsoid2(matter.updGround(), Transform(), body2, Transform());
    contacts.addBody(setIndex, ellipsoid1, ContactGeometry::Ellipsoid(radii1), center1);
    contacts.addBody(setIndex, ellipsoid2, ContactGeometry::Ellipsoid(radii2), center2);
    State state = system.realizeTopology();

    // Test a variety of positions.

    ellipsoid1.setQToFitTransform(state, Transform(Rotation(), Vec3(0))); // [-0.8, 0.8], [-1.7, 1.3], [-1.6, 2.6]
    ellipsoid2.setQToFitTransform(state, Transform(Rotation(), Vec3(2, 0, 0))); // [1.1, 3.1], [-1.2, 1.2], [-1.1, 1.7]
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    ellipsoid1.setQToFitTransform(state, Transform(Rotation(), Vec3(0))); // [-0.8, 0.8], [-1.7, 1.3], [-1.6, 2.6]
    ellipsoid2.setQToFitTransform(state, Transform(Rotation(), Vec3(1.5, 0, 0))); // [0.6, 2.6], [-1.2, 1.2], [-1.1, 1.7]
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 1);
    ASSERT(verifyEllipsoidContact(contacts.getContacts(state, setIndex)[0], radii1, radii2, center1, center2, ellipsoid1.getBodyTransform(state), ellipsoid2.getBodyTransform(state)));

    // Create a cloud of ellipsoids and find all contacts between them.

    MultibodySystem system2;
    SimbodyMatterSubsystem matter2(system2);
    GeneralContactSubsystem contacts2(system2);
    ContactSetIndex setIndex2 = contacts2.createContactSet();
    const int numEllipsoids = 100;
    for (int i = 0; i < numEllipsoids; i++) {
        MobilizedBody::Free ellipsoid(matter2.updGround(), Transform(), body, Transform());
        contacts2.addBody(setIndex2, ellipsoid, ContactGeometry::Ellipsoid(Vec3(0.1+random.getValue(), 0.1+random.getValue(), 0.1+random.getValue())), Vec3(0));
    }
    State state2 = system2.realizeTopology();
    for (MobilizedBodyIndex i(1); i <= numEllipsoids; i++) {
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(Vec3(random.getValue()*SimTK_PI, random.getValue()*SimTK_PI, random.getValue()*SimTK_PI));
        Vec3 pos = 5*Vec3(random.getValue(), random.getValue(), random.getValue());
        matter2.getMobilizedBody(i).setQToFitTransform(state2, Transform(rot, pos));
    }
    system2.realize(state2, Stage::Dynamics);

    // Verify each contact that was found.

    const Array_<Contact>& contact = contacts2.getContacts(state2, setIndex2);
    int errorCount = 0;
    for (int i = 0; i < (int) contact.size(); i++) {
        ASSERT(PointContact::isInstance(contact[i]));
        const PointContact& c = static_cast<const PointContact&>(contact[i]);
        const ContactGeometry::Ellipsoid& ellipsoid1 = reinterpret_cast<const ContactGeometry::Ellipsoid&>(contacts2.getBodyGeometry(setIndex2, c.getSurface1()));
        const ContactGeometry::Ellipsoid& ellipsoid2 = reinterpret_cast<const ContactGeometry::Ellipsoid&>(contacts2.getBodyGeometry(setIndex2, c.getSurface2()));
        const MobilizedBody& body1 = contacts2.getBody(setIndex2, c.getSurface1());
        const MobilizedBody& body2 = contacts2.getBody(setIndex2, c.getSurface2());
        if (!verifyEllipsoidContact(c, ellipsoid1.getRadii(), ellipsoid2.getRadii(), Vec3(0), Vec3(0), body1.getBodyTransform(state2), body2.getBodyTransform(state2)))
            errorCount++;
    }
    ASSERT(errorCount < (int)contact.size()/10);
}

/**
 * Check the set of faces in a contact.
 */

void verifyContactFaces(int* expected, int numExpected, const set<int>& found) {
    ASSERT(numExpected == found.size());
    for (int i = 0; i < numExpected; i++) {
        ASSERT(found.find(expected[i]) != found.end());
    }
}

void testHalfSpaceTriangleMesh() {
    // Create a triangle mesh consisting of two pyramids: one right side up and one upside down.
    
    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(0, 0, 1));
    vertices.push_back(Vec3(1, 0, 1));
    vertices.push_back(Vec3(0.5, 1, 0.5));
    vertices.push_back(Vec3(2, 1, 0));
    vertices.push_back(Vec3(2, 1, 1));
    vertices.push_back(Vec3(3, 1, 1));
    vertices.push_back(Vec3(3, 1, 0));
    vertices.push_back(Vec3(2.5, 0, 0.5));
    vector<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {0, 3, 4}, {3, 2, 4}, {2, 1, 4}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]+5);
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free b(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, b, mesh, Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::HalfSpace(), Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0, 1, 0))); // y < 1
    State state = system.realizeTopology();
    int bottomFaces[10] = {0, 1, 2, 3, 4, 5, 8, 9, 10, 11};
    for (Real depth = -0.25; depth < 1; depth += 0.1) {
        Vec3 center(0.1, 1-depth, 2.0);
        b.setQToFitTranslation(state, center);
        system.realize(state, Stage::Dynamics);
        const Array_<Contact>& contact = contacts.getContacts(state, setIndex);
        if (depth < 0.0) {
            ASSERT(contact.size() == 0);
        }
        else {
            ASSERT(contact.size() == 1);
            ASSERT(TriangleMeshContact::isInstance(contact[0]));
            const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contact[0]);
            ASSERT(c.getSurface1Faces().size() == 0);
            verifyContactFaces(bottomFaces, 10, c.getSurface2Faces());
         }
    }
}

void testSphereTriangleMesh() {
    // Create a triangle mesh consisting of two pyramids: one right side up and one upside down.
    
    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(0, 0, 1));
    vertices.push_back(Vec3(1, 0, 1));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(0.5, 1, 0.5));
    vertices.push_back(Vec3(2, 1, 0));
    vertices.push_back(Vec3(2, 1, 1));
    vertices.push_back(Vec3(3, 1, 1));
    vertices.push_back(Vec3(3, 1, 0));
    vertices.push_back(Vec3(2.5, 0, 0.5));
    vector<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {0, 3, 4}, {3, 2, 4}, {2, 1, 4}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]+5);
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free b(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, b, mesh, Transform());
    contacts.addBody(setIndex, matter.updGround(), ContactGeometry::Sphere(0.5), Transform(Vec3(0, 1, 0)));
    State state = system.realizeTopology();
    
    // Try various positions and make sure the results are correct.
    
    b.setQToFitTranslation(state, Vec3(0, -2, 0));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    b.setQToFitTranslation(state, Vec3(0, 1.51, 0));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    {
        b.setQToFitTranslation(state, Vec3(-0.5, 1.49, -0.5));
        system.realize(state, Stage::Dynamics);
        ASSERT(contacts.getContacts(state, setIndex).size() == 1);
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        int faces[] = {0, 1};
        verifyContactFaces(faces, 2, c.getSurface2Faces());
    }
    b.setQToFitTranslation(state, Vec3(-0.5, -0.51, -0.5));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    {
        b.setQToFitTranslation(state, Vec3(-0.5, -0.49, -0.5));
        system.realize(state, Stage::Dynamics);
        ASSERT(contacts.getContacts(state, setIndex).size() == 1);
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        int faces[] = {2, 3, 4, 5};
        verifyContactFaces(faces, 4, c.getSurface2Faces());
    }
    b.setQToFitTranslation(state, Vec3(-2.5, 1.51, -0.5));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    {
        b.setQToFitTranslation(state, Vec3(-2.5, 1.49, -0.5));
        system.realize(state, Stage::Dynamics);
        ASSERT(contacts.getContacts(state, setIndex).size() == 1);
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contacts.getContacts(state, setIndex)[0]);
        int faces[] = {8, 9, 10, 11};
        verifyContactFaces(faces, 4, c.getSurface2Faces());
    }
}

void testTriangleMeshTriangleMesh() {
    // Create two triangle meshes, each consisting of a pyramid.
    
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
    ContactGeometry::TriangleMesh mesh1(vertices, faceIndices);
    ContactGeometry::TriangleMesh mesh2(vertices, faceIndices);

    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralContactSubsystem contacts(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    ContactSetIndex setIndex = contacts.createContactSet();
    MobilizedBody::Free b1(matter.updGround(), Transform(), body, Transform());
    MobilizedBody::Free b2(matter.updGround(), Transform(), body, Transform());
    contacts.addBody(setIndex, b1, mesh1, Transform());
    contacts.addBody(setIndex, b2, mesh2, Transform());
    State state = system.realizeTopology();
    
    // Try some configurations that should not intersect.
    
    b1.setQToFitTranslation(state, Vec3(0));
    b2.setQToFitTranslation(state, Vec3(2));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    b1.setQToFitTranslation(state, Vec3(0));
    b2.setQToFitTranslation(state, Vec3(1.01, 0, 0));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    b1.setQToFitTranslation(state, Vec3(0));
    b2.setQToFitTranslation(state, Vec3(0, 1.01, 0));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    b1.setQToFitTranslation(state, Vec3(0));
    b2.setQToFitTranslation(state, Vec3(0, -1.01, 0));
    system.realize(state, Stage::Dynamics);
    ASSERT(contacts.getContacts(state, setIndex).size() == 0);
    
    // Now try ones that should intersect.
    
    int baseFaces[2] = {0, 1};
    int pointFaces[4] = {2, 3, 4, 5};
    {
        b1.setQToFitTranslation(state, Vec3(0));
        b2.setQToFitTranslation(state, Vec3(0, -0.99, 0));
        system.realize(state, Stage::Dynamics);
        Array_<Contact> contact = contacts.getContacts(state, setIndex);;
        ASSERT(contact.size() == 1);
        ASSERT(TriangleMeshContact::isInstance(contact[0]));
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contact[0]);
        if (contact[0].getSurface1() == 0) {
            verifyContactFaces(baseFaces, 2, c.getSurface1Faces());
            verifyContactFaces(pointFaces, 4, c.getSurface2Faces());
        }
        else {
            verifyContactFaces(pointFaces, 4, c.getSurface1Faces());
            verifyContactFaces(baseFaces, 2, c.getSurface2Faces());
        }
    }
    {
        b1.setQToFitTranslation(state, Vec3(0, -0.5, 0));
        b2.setQToFitTranslation(state, Vec3(0, 0.49, 0));
        system.realize(state, Stage::Dynamics);
        Array_<Contact> contact = contacts.getContacts(state, setIndex);;
        ASSERT(contact.size() == 1);
        ASSERT(TriangleMeshContact::isInstance(contact[0]));
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contact[0]);
        if (contact[0].getSurface1() == 0) {
            verifyContactFaces(pointFaces, 4, c.getSurface1Faces());
            verifyContactFaces(baseFaces, 2, c.getSurface2Faces());
        }
        else {
            verifyContactFaces(baseFaces, 2, c.getSurface1Faces());
            verifyContactFaces(pointFaces, 4, c.getSurface2Faces());
        }
    }
    {
        b1.setQToFitTranslation(state, Vec3(0.1, -0.5, 0));
        b2.setQToFitTranslation(state, Vec3(0, 0.49, 0.1));
        system.realize(state, Stage::Dynamics);
        Array_<Contact> contact = contacts.getContacts(state, setIndex);;
        ASSERT(contact.size() == 1);
        ASSERT(TriangleMeshContact::isInstance(contact[0]));
    }
    {
        b1.setQToFitTransform(state, Transform(Rotation(-0.5*Pi, ZAxis), Vec3(0, 0.5, 0)));
        b2.setQToFitTransform(state, Transform(Rotation(0.5*Pi, ZAxis), Vec3(1.9, -0.5, 0)));
        system.realize(state, Stage::Dynamics);
        Array_<Contact> contact = contacts.getContacts(state, setIndex);;
        ASSERT(contact.size() == 1);
        ASSERT(TriangleMeshContact::isInstance(contact[0]));
        const TriangleMeshContact& c = static_cast<const TriangleMeshContact&>(contact[0]);
        if (contact[0].getSurface1() == 0) {
            verifyContactFaces(pointFaces, 4, c.getSurface1Faces());
            verifyContactFaces(pointFaces, 4, c.getSurface2Faces());
        }
        else {
            verifyContactFaces(pointFaces, 4, c.getSurface1Faces());
            verifyContactFaces(pointFaces, 4, c.getSurface2Faces());
        }
    }
}

int main() {
    try {
        testHalfSpaceSphere();
        testSphereSphere();
        testHalfSpaceEllipsoid();
        testEllipsoidEllipsoid();
        testHalfSpaceTriangleMesh();
        testSphereTriangleMesh();
        testTriangleMeshTriangleMesh();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
