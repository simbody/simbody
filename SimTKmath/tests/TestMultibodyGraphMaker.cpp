/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-17 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

#include "SimTKmath.h"

using namespace SimTK;

void testSerialChain() {
    MultibodyGraphMaker mbgraph;
    mbgraph.addJointType("pin",  1);
    double nonzeroMass = 1.8;
    mbgraph.addBody("ground", 0, false);
    mbgraph.addBody("body1", nonzeroMass, false);
    mbgraph.addBody("body2", nonzeroMass, false);
    mbgraph.addBody("body3", nonzeroMass, false);
    mbgraph.addBody("body4", nonzeroMass, false);
    mbgraph.addBody("body5", nonzeroMass, false);

    mbgraph.addJoint("pin1", "pin", "ground", "body1", false);
    mbgraph.addJoint("pin2", "pin", "body1", "body2", false);
    mbgraph.addJoint("pin3", "pin", "body2", "body3", false);
    mbgraph.addJoint("pin4", "pin", "body3", "body4", false);
    mbgraph.addJoint("pin5", "pin", "body4", "body5", false);

    mbgraph.generateGraph();

    // First, some basic tests.
    SimTK_TEST(mbgraph.getNumJointTypes() == 3); // include weld and free.
    SimTK_TEST(mbgraph.getJointTypeNum("weld") == 0);
    SimTK_TEST(mbgraph.getJointTypeNum("free") == 1);
    SimTK_TEST(mbgraph.getJointTypeNum("pin") == 2); // 3rd joint.
    SimTK_TEST(mbgraph.getGroundBodyName() == "ground");
    SimTK_TEST(mbgraph.getNumBodies() == 6);
    SimTK_TEST(mbgraph.getBodyNum("body3") == 3); // ground is 0th body.
    SimTK_TEST(mbgraph.getNumJoints() == 5);
    SimTK_TEST(mbgraph.getJointNum("pin4") == 3);
    SimTK_TEST(mbgraph.getNumLoopConstraints() == 0);

    SimTK_TEST(mbgraph.getNumMobilizers() == 5);
    for (int i = 0; i < mbgraph.getNumMobilizers(); ++i) {
        const auto& mobilizer = mbgraph.getMobilizer(i);
        SimTK_TEST(!mobilizer.isAddedBaseMobilizer());
        SimTK_TEST(!mobilizer.isSlaveMobilizer());
        SimTK_TEST(!mobilizer.isReversedFromJoint());
        SimTK_TEST(mobilizer.getJointTypeName() == "pin");
        SimTK_TEST(mobilizer.getNumFragments() == 1);

        SimTK_TEST(mobilizer.getLevel() == i + 1);
    }
}

// This test ensures that the graph maker produces the same topology for a
// 5-body chain of pin joints (ground-body1-body2-body3-body4-body5) whether or
// not body2 is massless.
// This test was suggested by Sherm here:
// https://github.com/simbody/simbody/pull/592
void testIntermediateMasslessBody() {
    MultibodyGraphMaker mbgraph;
    mbgraph.addJointType("pin",  1);
    double nonzeroMass = 1.8;
    mbgraph.addBody("ground", 0, false);
    mbgraph.addBody("body1", nonzeroMass, false);
    mbgraph.addBody("body2", 0, false); // <-- Massless. 
    mbgraph.addBody("body3", nonzeroMass, false);
    mbgraph.addBody("body4", nonzeroMass, false);
    mbgraph.addBody("body5", nonzeroMass, false);

    mbgraph.addJoint("pin1", "pin", "ground", "body1", false);
    mbgraph.addJoint("pin2", "pin", "body1", "body2", false);
    mbgraph.addJoint("pin3", "pin", "body2", "body3", false);
    mbgraph.addJoint("pin4", "pin", "body3", "body4", false);
    mbgraph.addJoint("pin5", "pin", "body4", "body5", false);

    mbgraph.generateGraph();

    SimTK_TEST(mbgraph.getNumMobilizers() == 5);
    for (int i = 0; i < mbgraph.getNumMobilizers(); ++i) {
        const auto& mobilizer = mbgraph.getMobilizer(i);
        SimTK_TEST(!mobilizer.isAddedBaseMobilizer());
        SimTK_TEST(!mobilizer.isSlaveMobilizer());
        SimTK_TEST(!mobilizer.isReversedFromJoint());
        SimTK_TEST(mobilizer.getJointTypeName() == "pin");
        SimTK_TEST(mobilizer.getNumFragments() == 1);

        // Here's the test that the resulting graph is as expected:
        SimTK_TEST(mobilizer.getLevel() == i + 1);
    }
}

// Terminal massless bodies are not allowed (unless welded to a massful body).
void testTerminalMasslessBody() {
    MultibodyGraphMaker mbgraph;
    mbgraph.addJointType("pin",  1);
    mbgraph.addBody("ground", 0, false);
    mbgraph.addBody("body1", 0, false); // massless
    mbgraph.addJoint("pin1", "pin", "ground", "body1", false);
    SimTK_TEST_MUST_THROW_EXC(mbgraph.generateGraph(), std::runtime_error);

    // If the terminal massless body is welded, then there's no issue.
    mbgraph.deleteJoint("pin1");
    mbgraph.addJoint("pin1", "weld", "ground", "body1", false);
    mbgraph.generateGraph();
}

// This is a basic test for how MultibodyGraphMaker handles a system with a
// loop (four-bar linkage). There are many other settings surrounding loops
// that could be tested.
void testLoop() {
    MultibodyGraphMaker mbgraph;
    double nonzeroMass = 1.8;
    mbgraph.addJointType("pin",  1);
    mbgraph.addBody("ground", 0, false);
    mbgraph.addBody("body1", nonzeroMass, false);
    mbgraph.addBody("body2", nonzeroMass, false);
    mbgraph.addBody("body3", nonzeroMass, false);
    mbgraph.addJoint("pin1", "pin", "ground", "body1", false);
    mbgraph.addJoint("pin2", "pin", "body1", "body2", false);
    mbgraph.addJoint("pin3", "pin", "body2", "body3", false);
    mbgraph.addJoint("pin4", "pin", "body3", "ground", false);
    mbgraph.generateGraph();
    mbgraph.dumpGraph(std::cout);

    SimTK_TEST(mbgraph.getNumJointTypes() == 3); // include weld and free.
    SimTK_TEST(mbgraph.getJointTypeNum("weld") == 0);
    SimTK_TEST(mbgraph.getJointTypeNum("free") == 1);
    SimTK_TEST(mbgraph.getJointTypeNum("pin") == 2); // 3rd joint.
    SimTK_TEST(mbgraph.getGroundBodyName() == "ground");
    SimTK_TEST(mbgraph.getNumBodies() == 5); // body3 is split into two bodies.
    SimTK_TEST(mbgraph.getBody(4).name == "#body3_slave_1");
    SimTK_TEST(mbgraph.getNumJoints() == 4);
    // Still no loop constraints; the loop is broken by splitting a body and
    // adding a weld constraint.
    SimTK_TEST(mbgraph.getNumLoopConstraints() == 0);
    SimTK_TEST(mbgraph.getNumMobilizers() == 4);

    // Resulting mobilizers (breadth-first):
    // 0: pin1 ground->body1
    // 1: pin4 ground->body3 (reversed)
    // 2: pin2 body1->body2
    // 3: pin3 body2->#body3_slave_1

    SimTK_TEST(mbgraph.getMobilizer(0).getLevel() == 1);
    SimTK_TEST(mbgraph.getMobilizer(1).getLevel() == 1);
    SimTK_TEST(mbgraph.getMobilizer(2).getLevel() == 2);
    SimTK_TEST(mbgraph.getMobilizer(3).getLevel() == 3);
    for (int i = 0; i < mbgraph.getNumMobilizers(); ++i) {
        const auto& mobilizer = mbgraph.getMobilizer(i);
        SimTK_TEST(!mobilizer.isAddedBaseMobilizer());
        // The 4th mobilizer is a slave mobilizer (child body is a slave).
        SimTK_TEST(mobilizer.isSlaveMobilizer() == (i == 3));
        SimTK_TEST(mobilizer.isReversedFromJoint() == (i == 1));
        SimTK_TEST(mobilizer.getJointTypeName() == "pin");
        // Body 3 is split, and it's the outboard body for mobilizers 1 and 3.
        SimTK_TEST(mobilizer.getNumFragments() == (i == 1 || i == 3 ? 2 : 1));
    }
}

int main() {
    SimTK_START_TEST("TestMultibodyGraphMaker");
        SimTK_SUBTEST(testSerialChain);
        SimTK_SUBTEST(testIntermediateMasslessBody);
        SimTK_SUBTEST(testTerminalMasslessBody);
        SimTK_SUBTEST(testLoop);
    SimTK_END_TEST();
}
