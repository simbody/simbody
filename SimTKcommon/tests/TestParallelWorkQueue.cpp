/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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

class SetFlagTask : public ParallelWorkQueue::Task {
public:
    SetFlagTask(Array_<int>& flags, int index) : flags(flags), index(index) {
    }
    void execute() override {
        ASSERT(!flags[index]);
        flags[index] = true;
    }
private:
    Array_<int>& flags;
    int index;
};

void testParallelExecution() {
    const int numFlags = 500;
    Array_<int> flags(numFlags, false);
    ParallelWorkQueue queue(10);
    for (int i = 0; i < numFlags-10; i++)
        queue.addTask(new SetFlagTask(flags, i));
    queue.flush();
    for (int i = 0; i < numFlags; i++)
        ASSERT(flags[i] == (i < numFlags-10));
}

int main() {
    try {
        testParallelExecution();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
