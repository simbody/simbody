/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-17 Stanford University and the Authors.        *
 * Authors: Antoine Falisse                                                   *
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
#include "SimTKcommon/Testing.h"
#include <adolc\adolc.h>

#include <iostream>
using std::cout;
using std::endl;
using std::cin;

using namespace SimTK;

// Test derivative of simple function with ADOLC
void testAdouble() {
    double xp[1];
    xp[0] = -2.3;

    trace_on(1);
    adouble* x = new adouble[1];
    adouble* y = new adouble[1];
    x[0] <<= xp[0];
    y[0] = 3 * pow(x[0],3) + cos(x[0]) + 1;
    double y0[1];
    y[0] >>= y0[0];
    trace_off();

    double** J;
    J = myalloc(1,1);
    jacobian(1, 1, 1, xp, J);
    SimTK_TEST(J[0][0] == 48.355705212176716);

    delete[] y;
    delete[] x;
}

int main() {
    SimTK_START_TEST("TestADOLCCommon");
        SimTK_SUBTEST(testAdouble);
    SimTK_END_TEST();
}
