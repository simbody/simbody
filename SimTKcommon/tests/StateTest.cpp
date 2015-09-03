/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;
using std::string;

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

// Advance State by one stage from stage-1 to stage.
void advanceStage(State& state, Stage stage) {
    SimTK_TEST(state.getSystemStage() == stage.prev());
    for (SubsystemIndex sx(0); sx <state.getNumSubsystems(); ++sx)
        state.advanceSubsystemToStage(sx, stage);
    state.advanceSystemToStage(stage);
}

void testCacheValidity() {
    const SubsystemIndex Sub0(0), Sub1(1);
    State s;
    s.setNumSubsystems(2);

    SimTK_TEST(s.getSystemStage() == Stage::Empty);

    //-------------------------
    // "realize" Topology stage

    // Allocate at Topology stage a Model stage-invalidating state variable.
    const DiscreteVariableIndex dvx1TopoModel = 
        s.allocateDiscreteVariable(Sub1, Stage::Model, new Value<Real>(2));

    SimTK_TEST(s.getDiscreteVarAllocationStage(Sub1,dvx1TopoModel)
               == Stage::Topology);
    SimTK_TEST(s.getDiscreteVarInvalidatesStage(Sub1,dvx1TopoModel)
               == Stage::Model);
    SimTK_TEST(Value<Real>::downcast(s.getDiscreteVariable(Sub1,dvx1TopoModel))
               == Real(2));

    // Allocate at Topology stage a cache entry that depends on Model stage 
    // and is guaranteed to be valid at Time stage. In between (at Model or 
    // Instance stage) it *may* be valid if explicitly marked so.
    const CacheEntryIndex cx0TopoModel = s.allocateCacheEntry(Sub0, 
        Stage::Model, Stage::Time, new Value<int>(41));
    SimTK_TEST(s.getCacheEntryAllocationStage(Sub0,cx0TopoModel)
               == Stage::Topology);

    // Here is a cache entry allocated at Topology stage, with depends-on
    // Velocity and no good-by guarantee.
    const CacheEntryIndex cx0TopoVelocity = s.allocateCacheEntry(Sub0, 
        Stage::Velocity, Stage::Infinity, new Value<char>('v'));

    advanceStage(s, Stage::Topology);

    // Topology stage is realized.
    //----------------------------

    // Shouldn't be able to access cache entry here because this is less than
    // its "depends on" stage.
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx0TopoModel));

    //-------------------------
    // "realize" Model stage

    // Allocate at Model stage, a Position-invalidating state variable.
    const DiscreteVariableIndex dvx0ModelPos =
        s.allocateDiscreteVariable(Sub0, Stage::Position, new Value<int>(31));

    SimTK_TEST(s.getDiscreteVarAllocationStage(Sub0,dvx0ModelPos)
               == Stage::Model);
    SimTK_TEST(s.getDiscreteVarInvalidatesStage(Sub0,dvx0ModelPos)
               == Stage::Position);
    SimTK_TEST(Value<int>::downcast(s.getDiscreteVariable(Sub0,dvx0ModelPos))
               == 31);

    advanceStage(s, Stage::Model);
    // Model stage is realized.
    //----------------------------

    //----------------------------
    // "realize" Instance stage

    // Allocate a cache entry at Instance stage that has Time as depends-on
    // and also has a cross-subsystem dependency on discrete variable 
    // dvx0ModelPos.
    const CacheEntryIndex cx1InstanceTime = 
    s.allocateCacheEntryWithPrerequisites(Sub1, Stage::Time, Stage::Velocity,
        false, false, true, // depends on z
        {DiscreteVarKey(Sub0,dvx0ModelPos)}, 
        {CacheEntryKey(Sub0,cx0TopoModel)},
        new Value<string>("hasPrereqs_Time"));

    // Same but had depends-on Instance.
    const CacheEntryIndex cx1InstInst = 
    s.allocateCacheEntryWithPrerequisites(Sub1, Stage::Instance, Stage::Velocity,
        false, false, true, // depends on z
        {DiscreteVarKey(Sub0,dvx0ModelPos)}, 
        {CacheEntryKey(Sub0,cx0TopoModel)},
        new Value<string>("hasPrereqs_Instance"));

    // This attempt to create a cache-to-cache dependency should fail because
    // the prerequisite gets invalidated when Velocity stage changes but 
    // the "dependent" doesn't get invalidated unless Position stage does. Thus
    // a velocity could change, invlidating the prereq, but the downstream
    // cache entry still looks valid.
    SimTK_TEST_MUST_THROW(s.allocateCacheEntryWithPrerequisites(Sub1,
        Stage::Position, Stage::Infinity, false, false, false, {},
        {CacheEntryKey(Sub0,cx0TopoVelocity)}, new Value<int>(-1)));

    advanceStage(s, Stage::Instance);
    // Instance stage is realized.
    //----------------------------

    // Check that the dependency lists are right.
    CacheEntryKey ckey1(Sub1,cx1InstanceTime);
    CacheEntryKey ckey2(Sub1,cx1InstInst);
    SimTK_TEST(s.getZDependents().size() == 2);
    SimTK_TEST(s.getZDependents().contains(ckey1));
    SimTK_TEST(s.getZDependents().contains(ckey2));
    const DiscreteVarInfo& dvinfo = 
        s.getDiscreteVarInfo(DiscreteVarKey(Sub0,dvx0ModelPos));
    SimTK_TEST(dvinfo.getDependents().size() == 2);
    SimTK_TEST(dvinfo.getDependents().contains(ckey1));
    SimTK_TEST(dvinfo.getDependents().contains(ckey2));
    const CacheEntryInfo& ceinfo = 
        s.getCacheEntryInfo(CacheEntryKey(Sub0,cx0TopoModel));
    SimTK_TEST(ceinfo.getDependents().size() == 2);
    SimTK_TEST(ceinfo.getDependents().contains(ckey1));
    SimTK_TEST(ceinfo.getDependents().contains(ckey2));

    // Although cx0TopoModel *could* be valid at this point,
    // no one has said so, so we expect it to throw.
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx0TopoModel));

    // If we say it is valid, we should be able to obtain its value.
    s.markCacheValueRealized(Sub0, cx0TopoModel);
    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)) == 41);

    //----------------------------
    // "realize" Time stage
    advanceStage(s, Stage::Time);
    // Time stage is realized.
    //----------------------------

    // cx1InstanceTime isn't automatically valid but can be now.
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstanceTime));
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstInst));
    s.markCacheValueRealized(Sub1, cx1InstanceTime);
    s.markCacheValueRealized(Sub1, cx1InstInst);
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstanceTime)).get()
               == "hasPrereqs_Time");
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstInst)).get()
               == "hasPrereqs_Instance");

    // That cache entry does not depend on q or u so this should have no effect.
    s.updQ() = 0.;
    s.updU() = 0.;
    SimTK_TEST(s.getSystemStage() == Stage::Time); // unchanged
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstanceTime)).get()
               == "hasPrereqs_Time");

    // Changing prerequisites should make the cache entry
    // inaccessible again, although the stage should not change and the cache
    // entry should not get deallocated.
    s.updZ() = 0.; // z is prerequisite
    SimTK_TEST(s.getSystemStage() == Stage::Time); // unchanged
    SimTK_TEST(s.hasCacheEntry(CacheEntryKey(Sub1,cx1InstanceTime)));
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstanceTime));

    s.markCacheValueRealized(Sub1, cx1InstanceTime); // valid again
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstanceTime)).get()
               == "hasPrereqs_Time");

    // Modify design var prerequisite.
    Value<int>::updDowncast(s.updDiscreteVariable(Sub0,dvx0ModelPos)).upd()
        = 99;
    SimTK_TEST(s.hasCacheEntry(CacheEntryKey(Sub1,cx1InstanceTime)));
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstanceTime));
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstInst));
    SimTK_TEST(s.getSystemStage() == Stage::Time);

    s.markCacheValueRealized(Sub1, cx1InstanceTime); // valid again
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstanceTime)).get()
               == "hasPrereqs_Time");

    s.markCacheValueRealized(Sub1, cx1InstInst);
    SimTK_TEST(Value<string>::downcast(s.getCacheEntry(Sub1,cx1InstInst)).get()
               == "hasPrereqs_Instance");

    // Test state copying. Should copy through at least Instance stage.
    State s2(s); // copy construction

    SimTK_TEST(s2.getNumSubsystems() == s.getNumSubsystems());
    SimTK_TEST(s2.getSystemStage() >= Stage::Instance);
    SimTK_TEST(s2.hasCacheEntry(CacheEntryKey(Sub1,cx1InstanceTime)));
    SimTK_TEST(s2.hasCacheEntry(CacheEntryKey(Sub1,cx1InstInst)));

    // Dependency lists should have been reconstructed in the copy.
    SimTK_TEST(s2.getZDependents().size() == 2);
    SimTK_TEST(s2.getZDependents().contains(ckey1));
    SimTK_TEST(s2.getZDependents().contains(ckey2));
    const DiscreteVarInfo& dvinfo2 = 
        s2.getDiscreteVarInfo(DiscreteVarKey(Sub0,dvx0ModelPos));
    SimTK_TEST(dvinfo2.getDependents().size() == 2);
    SimTK_TEST(dvinfo2.getDependents().contains(ckey1));
    SimTK_TEST(dvinfo2.getDependents().contains(ckey2));
    const CacheEntryInfo& ceinfo2 = 
        s2.getCacheEntryInfo(CacheEntryKey(Sub0,cx0TopoModel));
    SimTK_TEST(ceinfo2.getDependents().size() == 2);
    SimTK_TEST(ceinfo2.getDependents().contains(ckey1));
    SimTK_TEST(ceinfo2.getDependents().contains(ckey2));


    s2.markCacheValueRealized(Sub1, cx1InstInst);
    SimTK_TEST(Value<string>::downcast(s2.getCacheEntry(Sub1,cx1InstInst)).get()
               == "hasPrereqs_Instance");

    // Invalidate cache entry prerequisite (modifying value is not enough).
    s.markCacheValueNotRealized(Sub0,cx0TopoModel);
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub1, cx1InstanceTime));
    SimTK_TEST(s.getSystemStage() == Stage::Time); // unchanged

    // Now modify a Model-stage state variable and realize again. This
    // should have invalidated Model stage and hence cx0TopoModel. It should
    // also have un-allocated the Instance-stage cache entry.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvx1TopoModel)) = 9;

    SimTK_TEST(s.getSystemStage() == Stage::Topology);
    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx0TopoModel));

    // Unallocating the cache entry should have removed it from its
    // prerequisite's dependency list.
    SimTK_TEST(!s.hasCacheEntry(CacheEntryKey(Sub1,cx1InstanceTime)));
    SimTK_TEST(s.getZDependents().empty());
    SimTK_TEST(ceinfo.getDependents().empty());

    //----------------------------
    // "realize" Model stage
    advanceStage(s, Stage::Model);
    // Model stage is realized.
    //----------------------------

    SimTK_TEST(!s.isCacheValueRealized(Sub0, cx0TopoModel));

    SimTK_TEST_MUST_THROW(s.getCacheEntry(Sub0, cx0TopoModel));

    // "calculate" the cache entry and mark it valid.
    Value<int>::updDowncast(s.updCacheEntry(Sub0,cx0TopoModel)) = 
        (int)(2*Value<Real>::downcast(s.getDiscreteVariable(Sub1,dvx1TopoModel)));
    s.markCacheValueRealized(Sub0, cx0TopoModel);

    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)) == 18);

    // Now modify the Model-stage variable again, but realize through
    // Time stage. We should be able to access the cache entry without
    // explicitly marking it valid.
    Value<Real>::updDowncast(s.updDiscreteVariable(Sub1, dvx1TopoModel)) = -100;
    advanceStage(s, Stage::Model);
    advanceStage(s, Stage::Instance);
    advanceStage(s, Stage::Time);

    SimTK_TEST(Value<int>::downcast(s.getCacheEntry(Sub0, cx0TopoModel)) == 18);

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
