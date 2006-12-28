/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

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
    s.advanceSubsystemToStage(0, Stage::Topology);
    s.advanceSystemToStage(Stage::Topology);

    Vector v3(3), v2(2);
    long q1 = s.allocateQ(0, v3);
    long q2 = s.allocateQ(0, v2);

    printf("q1,2=%d,%d\n", q1, q2);
    cout << s;

    long dv = s.allocateDiscreteVariable(0, Stage::Dynamics, new Value<int>(5));

    s.advanceSubsystemToStage(0, Stage::Model);
        //long dv2 = s.allocateDiscreteVariable(0, Stage::Position, new Value<int>(5));

    Value<int>::downcast(s.updDiscreteVariable(0, dv)) = 71;
    cout << s.getDiscreteVariable(0, dv) << endl;

    s.advanceSystemToStage(Stage::Model);

    cout << s;

  }
  catch(const std::exception& e) {
    printf("*** STATE TEST EXCEPTION\n%s\n***\n", e.what());
    return 1;
  }

    return 0;
}
