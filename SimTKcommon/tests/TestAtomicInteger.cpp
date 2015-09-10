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
    a = 3;
    a |= 4;
    ASSERT(a = 7);
    a &= 13;
    ASSERT(a == 5);
    a ^= 12;
    ASSERT(a == 9);
    a %= 7;
    ASSERT(a == 2);
    a <<= 2;
    ASSERT(a == 8);
    a >>= 1;
    ASSERT(a == 4);
}

void testParallelExecution() {
    AtomicInteger index;
    vector<int> flags(10000);
    for (int i = 0; i < (int)flags.size(); ++i)
        flags[i] = 0;
    ParallelExecutor executor;
    
    // See if the ++ operator is properly atomic.
    
    class SetFlagTask : public ParallelExecutor::Task {
    public:
        SetFlagTask(vector<int>& flags, AtomicInteger& index) : flags(flags), index(index) {
        }
        void execute(int i) override {
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
        for (int j = 0; j < (int)flags.size(); ++j)
            ASSERT(flags[j] == (j < 5000 ? i+1 : 0));
    }
    
    // See if the += operator is properly atomic.
    
    class IncrementTask : public ParallelExecutor::Task {
    public:
        IncrementTask(AtomicInteger& index) : index(index) {
        }
        void execute(int i) override {
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
    
    // See if the *= operator is properly atomic.
    
    class MultiplyTask : public ParallelExecutor::Task {
    public:
        MultiplyTask(AtomicInteger& index) : index(index) {
        }
        void execute(int i) override {
            index *= (i%500 == 0 ? 2 : -1);
        }
    private:
        AtomicInteger& index;
    };
    for (int i = 0; i < 100; ++i) {
        MultiplyTask task(index);
        index = 1;
        executor.execute(task, 4999);
        ASSERT(index == -1024);
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
