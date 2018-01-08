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

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

class SetFlagTask : public Parallel2DExecutor::Task {
public:
    SetFlagTask(Array_<Array_<int> >& flags, int& count) : flags(flags), count(count) {
    }
    void execute(int i, int j) override {
        flags[i][j]++;
        localCount++;
    }
    void initialize() override {
        localCount = 0;
    }
    void finish() override {
        count += localCount;
    }
private:
    Array_<Array_<int> >& flags;
    int& count;
    static thread_local int localCount;
};

/*static*/ thread_local int SetFlagTask::localCount = 0;

void clearFlags(Array_<Array_<int> >& flags) {
    int numFlags = flags.size();
    for (int i = 0; i < numFlags; ++i) {
        flags[i].resize(numFlags);
        for (int j = 0; j < numFlags; ++j)
            flags[i][j] = 0;
    }
}

void testParallelExecution() {
    int numFlags = 100;
    Parallel2DExecutor executor(numFlags);
    Array_<Array_<int> > flags(numFlags);
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        int count = 0;
        SetFlagTask task(flags, count);
        executor.execute(task, Parallel2DExecutor::FullMatrix);
        ASSERT(count == numFlags*numFlags);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == iter+1);
    }
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        int count = 0;
        SetFlagTask task(flags, count);
        executor.execute(task, Parallel2DExecutor::HalfMatrix);
        ASSERT(count == numFlags*(numFlags-1)/2);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == (j < i ? iter+1 : 0));
    }
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        int count = 0;
        SetFlagTask task(flags, count);
        executor.execute(task, Parallel2DExecutor::HalfPlusDiagonal);
        ASSERT(count == numFlags*(numFlags+1)/2);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == (j <= i ? iter+1 : 0));
    }
}

void testSingleThreadedExecution() {
    int numFlags = 100;
    Parallel2DExecutor executor(numFlags, 1);
    Array_<Array_<int> > flags(numFlags);
    int count = 0;
    SetFlagTask task(flags, count);
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::FullMatrix);
    ASSERT(count == numFlags*numFlags);
    for (int i = 0; i < numFlags; ++i)
        for (int j = 0; j < numFlags; ++j)
            ASSERT(flags[i][j] == 1);
    count = 0;
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::HalfMatrix);
    ASSERT(count == numFlags*(numFlags-1)/2);
    for (int i = 0; i < numFlags; ++i)
        for (int j = 0; j < numFlags; ++j)
            ASSERT(flags[i][j] == (j < i ? 1 : 0));
    count = 0;
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::HalfPlusDiagonal);
    ASSERT(count == numFlags*(numFlags+1)/2);
    for (int i = 0; i < numFlags; ++i)
        for (int j = 0; j < numFlags; ++j)
            ASSERT(flags[i][j] == (j <= i ? 1 : 0));
}

int main() {
    try {
        testParallelExecution();
        testSingleThreadedExecution();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
