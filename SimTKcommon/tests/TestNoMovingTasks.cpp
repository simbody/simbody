/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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

#include "SimTKcommon.h"

#include <iostream>
#include <thread>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using namespace SimTK;

static int num_cpus;

class TestTask : public Parallel2DExecutor::Task {
public:
    TestTask() {}
    void execute(int, int) override {
        ASSERT(threadId == std::this_thread::get_id());
    }
    void initialize() override {
        threadId = std::this_thread::get_id();
    }
    void finish() override {
        ASSERT(threadId == std::this_thread::get_id());

        threadId = std::thread::id{};
    }

private:
    static thread_local std::thread::id threadId;
    static thread_local int counter;
};

thread_local std::thread::id TestTask::threadId;

void testDefault() {
    for (int size = 1; size < num_cpus * 2 + 1; size++) {
        Parallel2DExecutor executor(size);
        TestTask task;
        executor.execute(task, Parallel2DExecutor::FullMatrix);
    }
}

void testFixedMultiple()
{
    for (int size = 1; size < num_cpus * 2 + 1; size++) {
        ParallelExecutor e(num_cpus > 1 ? num_cpus : 1);
        Parallel2DExecutor executor(size, e);
        TestTask task;
        executor.execute(task, Parallel2DExecutor::FullMatrix);
    }
}

void testFixedSingle()
{
    for (int size = 1; size < num_cpus * 2 + 1; size++) {
        ParallelExecutor e(1);
        Parallel2DExecutor executor(size, e);
        TestTask task;
        executor.execute(task, Parallel2DExecutor::FullMatrix);
    }
}

int main() {
    SimTK_START_TEST("TestStableThreadId");
        num_cpus = ParallelExecutor::getNumProcessors();
        SimTK_SUBTEST(testDefault);
        SimTK_SUBTEST(testFixedMultiple);
        SimTK_SUBTEST(testFixedSingle);
    SimTK_END_TEST();
    return 0;
}
