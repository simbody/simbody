/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

using namespace SimTK;

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


  try {
    State s;
    s.setNSubsystems(1);
    s.advanceSubsystemToStage(SubsystemId(0), Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    Vector v3(3), v2(2);
    int q1 = s.allocateQ(SubsystemId(0), v3);
    int q2 = s.allocateQ(SubsystemId(0), v2);

    int e1 = s.allocateEvent(SubsystemId(0), Stage::Position, 3);
    int e2 = s.allocateEvent(SubsystemId(0), Stage::Instance, 2);

    printf("q1,2=%d,%d\n", q1, q2);
    printf("e1,2=%d,%d\n", e1, e2);

    cout << s;

    long dv = s.allocateDiscreteVariable(SubsystemId(0), Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(SubsystemId(0), Stage::Model);
        //long dv2 = s.allocateDiscreteVariable(SubsystemId(0), Stage::Position, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(SubsystemId(0), dv)) = 71;
    cout << s.getDiscreteVariable(SubsystemId(0), dv) << endl;

    // Can't ask for the time before Stage::Model, but if you could it would be NaN.
    // Advancing to Stage::Model sets t=0.
    //cout << "BEFORE ADVANCE TO MODEL, t=" << s.getTime() << endl;

    s.advanceSystemToStage(Stage::Model);

    cout << "AFTER ADVANCE TO MODEL, t=" << s.getTime() << endl;

    printf("nevents=%d, by stage:\n", s.getNEvents());
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage::getValue(j);
        cout << g.getName() << ": " << s.getNEventsByStage(g) << endl;
    }

    printf("subsys 0 by stage:\n");
    for (int j=0; j<Stage::NValid; ++j) {
        Stage g = Stage::getValue(j);
        cout << g.getName() << ": " << s.getNEventsByStage(SubsystemId(0),g) << endl;
    }
    cout << "State s=" << s;

    s.clear();
    cout << "after clear(), State s=" << s;

  }
  catch(const std::exception& e) {
    printf("*** STATE TEST EXCEPTION\n%s\n***\n", e.what());
    return 1;
  }

    return 0;
}
