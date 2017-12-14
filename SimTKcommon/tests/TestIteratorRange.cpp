/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include <map>

using namespace SimTK;

void testFirstLast() {
    std::vector<int> v {1, 5, 6, 7, 9, 10, 15};
    
    {
        std::vector<int> indices {6, 7, 9};
        int index = 0;
        for (auto& x : makeIteratorRange(
                    std::lower_bound(v.begin(), v.end(), 6),
                    std::upper_bound(v.begin(), v.end(), 9))) {
            SimTK_TEST(indices[index] == x);
            index++;
        }
    }
}

void testPair() {
    std::multimap<std::string, int> map;
    map.insert({"avocado", 5});
    map.insert({"broccoli", 5});
    map.insert({"broccoli", 16});
    map.insert({"eggplant", 1});
    map.insert({"guineafowl", 11});
    map.insert({"broccoli", 17});
    map.insert({"broccoli", 5});
    map.insert({"avocado", 71});
    map.insert({"jackfruit", 8});

    {
        std::vector<int> avocado {5, 71};
        int index = 0;
        for (auto& x : makeIteratorRange(map.equal_range("avocado"))) {
            SimTK_TEST(avocado[index] == x.second);
            index++;
        }
    }

    {
        std::vector<int> broccoli {5, 16, 17, 5};
        int index = 0;
        for (auto& x : makeIteratorRange(map.equal_range("broccoli"))) {
            SimTK_TEST(broccoli[index] == x.second);
            index++;
        }
    }
}

int main() {
    SimTK_START_TEST("TestIteratorRange");

        SimTK_SUBTEST(testFirstLast);
        SimTK_SUBTEST(testPair);

    SimTK_END_TEST();
}
