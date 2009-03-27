/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * Run the State class few some paces.
 */

#include "SimTKcommon.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS((cond), "Assertion failed.");}


void testLowestModified() {
    const SubsystemIndex Sub0(0), Sub1(1);
    State s;
    s.setNSubsystems(2);
    ASSERT(s.getSystemStage()==Stage::Empty);
    ASSERT(s.getSubsystemStage(Sub0)==Stage::Empty && s.getSubsystemStage(Sub1)==Stage::Empty);
    ASSERT(s.getLowestStageModified()==Stage::Topology);

    const DiscreteVariableIndex dvxModel = s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // "realize" Topology stage
    s.advanceSubsystemToStage(Sub0, Stage::Topology);
    s.advanceSubsystemToStage(Sub1, Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);
    ASSERT(s.getLowestStageModified()==Stage::Topology);    // shouldn't have changed

    const DiscreteVariableIndex dvxInstance = s.allocateDiscreteVariable(Sub0, Stage::Instance, new Value<int>(-4));

    // A Model-stage variable must be allocated before Topology is realized, and this condition
    // should be tested even in Release mode.
    try {
        s.allocateDiscreteVariable(Sub0, Stage::Model, new Value<int>(0));
        ASSERT(!"Shouldn't have allowed allocation of Model-stage var here");
    } catch (...) {}

    // "realize" Model stage
    s.advanceSubsystemToStage(Sub0, Stage::Model);
    s.advanceSubsystemToStage(Sub1, Stage::Model);
    s.advanceSystemToStage(Stage::Model);
    ASSERT(s.getSystemStage() == Stage::Model);

    ASSERT(s.getLowestStageModified()==Stage::Topology); // shouldn't have changed
    s.resetLowestStageModified();
    ASSERT(s.getLowestStageModified()==Stage::Instance);  // i.e., lowest invalid stage

    // This variable invalidates Instance stage, so shouldn't change anything now.
    ASSERT(Value<int>::downcast(s.getDiscreteVariable(Sub0, dvxInstance))==-4);
    Value<int>::downcast(s.updDiscreteVariable(Sub0, dvxInstance)) = 123;
    ASSERT(Value<int>::downcast(s.getDiscreteVariable(Sub0, dvxInstance))== 123);

    ASSERT(s.getSystemStage()==Stage::Model);
    ASSERT(s.getLowestStageModified()==Stage::Instance);

    // This variable invalidates Model Stage, so should back up the stage to Topology,
    // invalidating Model.
    Value<Real>::downcast(s.updDiscreteVariable(Sub1, dvxModel)) = 29;
    ASSERT(s.getSubsystemStage(Sub1)==Stage::Topology);
    ASSERT(s.getLowestStageModified()==Stage::Model);

    // Now realize Model stage again; shouldn't affect lowestStageModified.
    s.advanceSubsystemToStage(Sub0, Stage::Model);
    s.advanceSubsystemToStage(Sub1, Stage::Model);
    s.advanceSystemToStage(Stage::Model);
    ASSERT(s.getSystemStage() == Stage::Model);

    ASSERT(s.getLowestStageModified()==Stage::Model); // shouldn't have changed
}

int main() {
    int major,minor,build;
    char out[100];
    const char* keylist[] = { "version", "library", "type", "debug", "authors", "copyright", "svn_revision", 0 };

    SimTK_version_SimTKcommon(&major,&minor,&build);
    std::printf("==> SimTKcommon library version: %d.%d.%d\n", major, minor, build);
    std::printf("    SimTK_about_SimTKcommon():\n");
    for (const char** p = keylist; *p; ++p) {
        SimTK_about_SimTKcommon(*p, 100, out);
        std::printf("      about(%s)='%s'\n", *p, out);
    }


  try {
    testLowestModified();

    State s;
    s.setNSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);

    // Can't ask for the time before Stage::Topology, but if you could it would be NaN.
    s.advanceSystemToStage(Stage::Topology);

    // Advancing to Stage::Topology sets t=0.
    cout << "AFTER ADVANCE TO TOPOLOGY, t=" << s.getTime() << endl;

    ASSERT(s.getTime()==0);

    Vector v3(3), v2(2);
    QIndex q1 = s.allocateQ(SubsystemIndex(0), v3);
    QIndex q2 = s.allocateQ(SubsystemIndex(0), v2);

    EventTriggerByStageIndex e1 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Position, 3);
    EventTriggerByStageIndex e2 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Instance, 2);

    printf("q1,2=%d,%d\n", (int)q1, (int)q2);
    printf("e1,2=%d,%d\n", (int)e1, (int)e2);

    cout << s;

    DiscreteVariableIndex dv = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Model);
        //long dv2 = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Position, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(SubsystemIndex(0), dv)) = 71;
    cout << s.getDiscreteVariable(SubsystemIndex(0), dv) << endl;


    s.advanceSystemToStage(Stage::Model);

    cout << "AFTER ADVANCE TO MODEL, t=" << s.getTime() << endl;

    printf("ntriggers=%d, by stage:\n", s.getNEventTriggers());
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage::getValue(j);
        cout << g.getName() << ": " << s.getNEventTriggersByStage(g) << endl;
    }

    printf("subsys 0 by stage:\n");
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage::getValue(j);
        cout << g.getName() << ": " << s.getNEventTriggersByStage(SubsystemIndex(0),g) << endl;
    }
    cout << "State s=" << s;

    s.clear();
    cout << "after clear(), State s=" << s;

  }
  catch(const std::exception& e) {
    printf("*** STATE TEST EXCEPTION\n%s\n***\n", e.what());
    return 1;
  }

    return 0;
}
