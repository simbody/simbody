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
        
        const vector<Contact>& contact = contacts.getContacts(state, setIndex);
        if (centerInGround[1] > radius+1) {
            ASSERT(contact.size() == 0);
        }
        else {
            ASSERT(contact.size() == 1);
            assertEqual(contact[0].getFirstBody(), 1);
            assertEqual(contact[0].getSecondBody(), 0);
            assertEqual(contact[0].getNormal(), Vec3(0, 1, 0));
            Real depth = radius-centerInGround[1]+1;
            assertEqual(contact[0].getDepth(), depth);
            assertEqual(contact[0].getRadius(), std::sqrt(radius*depth));
            assertEqual(contact[0].getLocation(), Vec3(centerInGround[0], 1-0.5*depth, centerInGround[2]));
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
        
        const vector<Contact>& contact = contacts.getContacts(state, setIndex);
        for (int i = 0; i < contact.size(); i++) {
            int body1 = contact[i].getFirstBody();
            int body2 = contact[i].getSecondBody();
            Vec3 delta = centerInGround[body2]-centerInGround[body1];
            assertEqual(delta.normalize(), contact[i].getNormal());
            assertEqual(delta.norm(), radius[body1]+radius[body2]-contact[i].getDepth());
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

int main() {
    try {
        testHalfSpaceSphere();
        testSphereSphere();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
