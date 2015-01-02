/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * Run the State class few some paces.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;




//void testLowestModified() {
//    const SubsystemIndex Sub0(0), Sub1(1);
//    State s;
//    s.setNumSubsystems(2);
//    SimTK_TEST(s.getSystemStage()==Stage::Empty);
//    SimTK_TEST(s.getSubsystemStage(Sub0)==Stage::Empty && s.getSubsystemStage(Sub1)==Stage::Empty);
//    SimTK_TEST(s.getLowestStageModified()==Stage::Topology);
//
//    const DiscreteVariableIndex dvxModel = s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));
//
//    // "realize" Topology stage
//    s.advanceSubsystemToStage(Sub0, Stage::Topology);
//    s.advanceSubsystemToStage(Sub1, Stage::Topology);
//    s.advanceSystemToStage(Stage::Topology);
//    SimTK_TEST(s.getLowestStageModified()==Stage::Topology);    // shouldn't have changed
//
//    const DiscreteVariableIndex dvxInstance = s.allocateDiscreteVariable(Sub0, Stage::Instance, new Value<int>(-4));
//
//    // A Model-stage variable must be allocated before Topology is realized, and this condition
//    // should be tested even in Release mode.
//    try {
//        s.allocateDiscreteVariable(Sub0, Stage::Model, new Value<int>(0));
//        SimTK_TEST(!"Shouldn't have allowed allocation of Model-stage var here");
//    } catch (...) {}
//
//    // "realize" Model stage
//    s.advanceSubsystemToStage(Sub0, Stage::Model);
//    s.advanceSubsystemToStage(Sub1, Stage::Model);
//    s.advanceSystemToStage(Stage::Model);
//    SimTK_TEST(s.getSystemStage() == Stage::Model);
//
//    SimTK_TEST(s.getLowestStageModified()==Stage::Topology); // shouldn't have changed
//    s.resetLowestStageModified();
//    SimTK_TEST(s.getLowestStageModified()==Stage::Instance);  // i.e., lowest invalid stage
//
//    // This variable invalidates Instance stage, so shouldn't change anything now.
//    SimTK_TEST(Value<int>::downcast(s.getDiscreteVariable(Sub0, dvxInstance))==-4);
//    Value<int>::downcast(s.updDiscreteVariable(Sub0, dvxInstance)) = 123;
//    SimTK_TEST(Value<int>::downcast(s.getDiscreteVariable(Sub0, dvxInstance))== 123);
//
//    SimTK_TEST(s.getSystemStage()==Stage::Model);
//    SimTK_TEST(s.getLowestStageModified()==Stage::Instance);
//
//    // This variable invalidates Model Stage, so should back up the stage to Topology,
//    // invalidating Model.
//    Value<Real>::downcast(s.updDiscreteVariable(Sub1, dvxModel)) = 29;
//    SimTK_TEST(s.getSubsystemStage(Sub1)==Stage::Topology);
//    SimTK_TEST(s.getLowestStageModified()==Stage::Model);
//
//    // Now realize Model stage again; shouldn't affect lowestStageModified.
//    s.advanceSubsystemToStage(Sub0, Stage::Model);
//    s.advanceSubsystemToStage(Sub1, Stage::Model);
//    s.advanceSystemToStage(Stage::Model);
//    SimTK_TEST(s.getSystemStage() == Stage::Model);
//
//    SimTK_TEST(s.getLowestStageModified()==Stage::Model); // shouldn't have changed
//}

void testCacheValidity() {
    const SubsystemIndex Sub0(0), Sub1(1);
    State s;
    s.setNumSubsystems(2);

    // Allocate a Model stage-invalidating state variable.
    const DiscreteVariableIndex dvxModel = s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    // This cache entry depends on Model stage state and is guaranteed to be valid at Time stage.
    // In between (at Model or Instance stage) it *may* be valid if explicitly marked so.
    const CacheEntryIndex cx = s.allocateCacheEntry(Sub0, 
        Stage::Model, Stage::Time, new Value<int>(41));

    // "realize" Topology stage
    s.advanceSubsystemToStage(Sub0, Stage::Topology);
    s.advanceSubsystemToStage(Sub1, Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    // Shouldn't be able to access cache entry here because this is less than its
    // "depends on" stage.
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx));


    // "realize" Model stage
    s.advanceSubsystemToStage(Sub0, Stage::Model);
    s.advanceSubsystemToStage(Sub1, Stage::Model);
    s.advanceSystemToStage(Stage::Model);

    // "realize" Instance stage
    s.advanceSubsystemToStage(Sub0, Stage::Instance);
    s.advanceSubsystemToStage(Sub1, Stage::Instance);
    s.advanceSystemToStage(Stage::Instance);

    // Although the cache entry *could* be valid at this point,
    // no one has said so, so we expect it to throw.
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx));

    // If we say it is valid, we should be able to obtain its value.
    s.markCacheValueRealized(Sub0, cx);
    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx)) == 41);

    // Now modify a Model-stage state variable and realize again. This
    // shoudl have invalidated the cache entry.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvxModel)) = 9;
    // "realize" Model stage
    s.advanceSubsystemToStage(Sub0, Stage::Model);
    s.advanceSubsystemToStage(Sub1, Stage::Model);
    s.advanceSystemToStage(Stage::Model);

    SimTK_TEST(!s.isCacheValueRealized(Sub0, cx));

    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx));

    // "calculate" the cache entry and mark it valid.
    Value<int>::updDowncast(s.updCacheEntry(Sub0,cx)) = 
        (int)(2 * Value<Real>::downcast(s.getDiscreteVariable(Sub1, dvxModel)));
    s.markCacheValueRealized(Sub0, cx);

    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx)) == 18);

    // Now modify the Model-stage variable again, but realize through
    // Time stage. We should be able to access the cache entry without
    // explicitly marking it valid.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvxModel)) = -100;
    s.advanceSubsystemToStage(Sub0, Stage::Model);
    s.advanceSubsystemToStage(Sub1, Stage::Model);
    s.advanceSystemToStage(Stage::Model);
    s.advanceSubsystemToStage(Sub0, Stage::Instance);
    s.advanceSubsystemToStage(Sub1, Stage::Instance);
    s.advanceSystemToStage(Stage::Instance);
    s.advanceSubsystemToStage(Sub0, Stage::Time);
    s.advanceSubsystemToStage(Sub1, Stage::Time);
    s.advanceSystemToStage(Stage::Time);

    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx)) == 18);

}

void testMisc() {
    State s;
    s.setNumSubsystems(1);
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Topology);

    // Can't ask for the time before Stage::Topology, but if you could it would be NaN.
    s.advanceSystemToStage(Stage::Topology);

    // Advancing to Stage::Topology sets t=0.
    cout << "AFTER ADVANCE TO TOPOLOGY, t=" << s.getTime() << endl;

    SimTK_TEST(s.getTime()==0);

    Vector v3(3), v2(2);
    QIndex q1 = s.allocateQ(SubsystemIndex(0), v3);
    QIndex q2 = s.allocateQ(SubsystemIndex(0), v2);

    EventTriggerByStageIndex e1 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Position, 3);
    EventTriggerByStageIndex e2 = s.allocateEventTrigger(SubsystemIndex(0), Stage::Instance, 2);

    printf("q1,2=%d,%d\n", (int)q1, (int)q2);
    printf("e1,2=%d,%d\n", (int)e1, (int)e2);

    //cout << s;

    DiscreteVariableIndex dv = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Model);
        //long dv2 = s.allocateDiscreteVariable(SubsystemIndex(0), Stage::Position, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(SubsystemIndex(0), dv)) = 71;
    cout << s.getDiscreteVariable(SubsystemIndex(0), dv) << endl;


    s.advanceSystemToStage(Stage::Model);

    cout << "AFTER ADVANCE TO MODEL, t=" << s.getTime() << endl;

    // Event triggers are available at Instance stage.
    s.advanceSubsystemToStage(SubsystemIndex(0), Stage::Instance);
    s.advanceSystemToStage(Stage::Instance);

    printf("ntriggers=%d, by stage:\n", s.getNEventTriggers());
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage(j);
        cout << g.getName() << ": " << s.getNEventTriggersByStage(g) << endl;
    }

    printf("subsys 0 by stage:\n");
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage(j);
        cout << g.getName() << ": " << s.getNEventTriggersByStage(SubsystemIndex(0),g) << endl;
    }
    //cout << "State s=" << s;

    s.clear();
    //cout << "after clear(), State s=" << s;
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


    SimTK_START_TEST("StateTest");
        //SimTK_SUBTEST(testLowestModified);
        SimTK_SUBTEST(testCacheValidity);
        SimTK_SUBTEST(testMisc);
    SimTK_END_TEST();
}
