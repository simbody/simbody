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

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

void testOperators() {
    AtomicInteger a = 5;
    ASSERT(a == 5);
    ASSERT(a != 6);
    a = 6;
    ASSERT(a == 6);
    ASSERT(++a == 7);
    ASSERT(a == 7);
    ASSERT(a++ == 7);
    ASSERT(a == 8);
    ASSERT(--a == 7);
    ASSERT(a == 7);
    ASSERT(a-- == 7);
    ASSERT(a == 6);
    a += 2;
    ASSERT(a == 8);
    a -= 5;
    ASSERT(a == 3);
    a *= -5;
    ASSERT(a == -15);
    a /= 3;
    ASSERT(a == -5);
    int i = a;
    ASSERT(i == -5);
}

void testParallelExecution() {
    AtomicInteger index;
    vector<int> flags(10000);
    for (int i = 0; i < flags.size(); ++i)
        flags[i] = 0;
    ParallelExecutor executor;
    
    // See if the ++ operator is properly atomic.
    
    class SetFlagTask : public ParallelExecutor::Task {
    public:
        SetFlagTask(vector<int>& flags, AtomicInteger& index) : flags(flags), index(index) {
        }
        void execute(int i) {
            flags[index++]++;
        }
    private:
        vector<int>& flags;
        AtomicInteger& index;
    };
    for (int i = 0; i < 100; ++i) {
        SetFlagTask task(flags, index);
        index = 0;
        executor.execute(task, 5000);
        ASSERT(index == 5000);
        for (int j = 0; j < flags.size(); ++j)
            ASSERT(flags[j] == (j < 5000 ? i+1 : 0));
    }
    
    // See if the += operator is properly atomic.
    
    class IncrementTask : public ParallelExecutor::Task {
    public:
        IncrementTask(AtomicInteger& index) : index(index) {
        }
        void execute(int i) {
            index += 2;
        }
    private:
        AtomicInteger& index;
    };
    for (int i = 0; i < 100; ++i) {
        IncrementTask task(index);
        index = 0;
        executor.execute(task, 5000);
        ASSERT(index == 10000);
    }
}

int main() {
    try {
        testOperators();
        testParallelExecution();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
