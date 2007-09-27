/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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
using std::sqrt;
using std::vector;
using namespace SimTK;

class Color : public TypesafeEnum<Color> {
public:
    static const Color Red;
    static const Color Green;
    static const Color Blue;
private:
    Color() : TypesafeEnum<Color>() {
    }
};

const Color Color::Red;
const Color Color::Green;
const Color Color::Blue;

int main() {
    try {
        const vector<TypesafeEnum<Color> > values = Color::getAllValues();
        assert(values.size() == 3);
        assert(values[Color::Red.getIndex()] == Color::Red);
        EnumSet<Color> set1;
        assert(set1.size() == 0);
        set1 += Color::Red;
        set1 += Color::Blue;
        assert(set1.size() == 2);
        assert(set1.contains(Color::Red));
        assert(set1.contains(Color::Blue));
        assert(!set1.contains(Color::Green));
        EnumSet<Color> set2;
        assert(set2.size() == 0);
        set2 += Color::Blue;
        assert(set2.size() == 1);
        assert(!set2.contains(Color::Red));
        assert(set2.contains(Color::Blue));
        assert(!set2.contains(Color::Green));
        assert(set1.containsAll(set2));
        assert(!set2.containsAll(set1));
        assert(set1.containsAny(set2));
        assert(set2.containsAny(set1));
        assert(set1 != set2);
        set1 -= Color::Red;
        assert(set1.size() == 1);
        assert(set1 == set2);
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
