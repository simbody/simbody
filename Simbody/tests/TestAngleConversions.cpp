/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

/**
 * Test converting a SimbodyMatterSubsystem between quaternion and Euler angle representations.
 */

int main() {

    // Build a system consisting of a chain of bodies with every possible mobilizer.

    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)));
    body.addDecoration(DecorativeSphere(.1));
    Random::Uniform random(0.0, 2.0);
    MobilizedBody lastBody = MobilizedBody::Pin(matter.Ground(), Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Slider(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Universal(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Cylinder(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::BendStretch(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Planar(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Gimbal(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Ball(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Translation(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Free(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::LineOrientation(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::FreeLine(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Weld(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Screw(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())), 0.5);
    lastBody = MobilizedBody::Ellipsoid(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    mbs.realizeTopology();
    State& s = mbs.updDefaultState();
    mbs.realizeModel(s);

    // Choose a random initial conformation.

    for (int i = 0; i < s.getNQ(); ++i)
        s.updQ()[i] = random.getValue();
    mbs.realize(s, Stage::Instance);
    // The only constraints are the quaternions -- normalize them.
    Vector temp;
    mbs.project(s, 0.01);
    mbs.realize(s, Stage::Position);

    // Convert to Euler angles and make sure the positions are all the same.

    State euler = s;
    matter.convertToEulerAngles(s, euler);
    mbs.realize(euler, Stage::Position);
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (body.getBodyOriginLocation(euler)-body.getBodyOriginLocation(s)).norm();
        ASSERT(dist < 1e-5);
    }

    // Now convert back to quaternions and make sure the positions are still the same.

    State quaternions = s;
    matter.convertToQuaternions(euler, quaternions);
    mbs.realize(quaternions, Stage::Position);
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (body.getBodyOriginLocation(quaternions)-body.getBodyOriginLocation(s)).norm();
        ASSERT(dist < 1e-5);
    }

    // Compare the state variables to see if they have been accurately reproduced.

    mbs.project(s, 0.01); // Normalize the quaternions
    Real diff = std::sqrt((s.getQ()-quaternions.getQ()).normSqr()/s.getNQ());
    ASSERT(diff < 1e-5);

    std::cout << "Done" << std::endl;
}
