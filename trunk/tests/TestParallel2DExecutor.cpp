/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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

#include "SimTKcommon.h"

#include <iostream>
#include <vector>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

class SetFlagTask : public Parallel2DExecutor::Task {
public:
    SetFlagTask(vector<vector<int> >& flags) : flags(flags) {
    }
    void execute(int i, int j) {
        flags[i][j]++;
    }
private:
    vector<vector<int> >& flags;
};

void clearFlags(vector<vector<int> >& flags) {
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
    vector<vector<int> > flags(numFlags);
    SetFlagTask task(flags);
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        executor.execute(task, Parallel2DExecutor::FullMatrix);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == iter+1);
    }
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        executor.execute(task, Parallel2DExecutor::HalfMatrix);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == (j < i ? iter+1 : 0));
    }
    clearFlags(flags);
    for (int iter = 0; iter < 100; ++iter) {
        executor.execute(task, Parallel2DExecutor::HalfPlusDiagonal);
        for (int i = 0; i < numFlags; ++i)
            for (int j = 0; j < numFlags; ++j)
                ASSERT(flags[i][j] == (j <= i ? iter+1 : 0));
    }
}

void testSingleThreadedExecution() {
    int numFlags = 100;
    Parallel2DExecutor executor(numFlags, 1);
    vector<vector<int> > flags(numFlags);
    SetFlagTask task(flags);
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::FullMatrix);
    for (int i = 0; i < numFlags; ++i)
        for (int j = 0; j < numFlags; ++j)
            ASSERT(flags[i][j] == 1);
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::HalfMatrix);
    for (int i = 0; i < numFlags; ++i)
        for (int j = 0; j < numFlags; ++j)
            ASSERT(flags[i][j] == (j < i ? 1 : 0));
    clearFlags(flags);
    executor.execute(task, Parallel2DExecutor::HalfPlusDiagonal);
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
