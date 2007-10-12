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
    Color(int index, char* name) : TypesafeEnum<Color>(index, name) {
    }
    static void initValues() {
        new(&const_cast<Color&>(Red)) Color(0, "Red");
        new(&const_cast<Color&>(Green)) Color(1, "Green");
        new(&const_cast<Color&>(Blue)) Color(2, "Blue");
    }
    friend class TypesafeEnum<Color>;
};

Color acolor = Color::Green; // Verify that accessing a constant that has not yet been initialized works correctly.

const Color Color::Red;
const Color Color::Green;
const Color Color::Blue;

void verifyContents(EnumSet<Color> set, bool hasRed, bool hasBlue, bool hasGreen) {
    ASSERT(set.contains(Color::Red) == hasRed);
    ASSERT(set.contains(Color::Blue) == hasBlue);
    ASSERT(set.contains(Color::Green) == hasGreen);
    int size = 0;
    if (hasRed)
        size++;
    if (hasBlue)
        size++;
    if (hasGreen)
        size++;
    ASSERT(set.size() == size);
    bool redFound = false, blueFound = false, greenFound = false;
    for (EnumSet<Color>::iterator iter = set.begin(); iter != set.end(); ++iter) {
        if (*iter == Color::Red) {
            ASSERT(!redFound);
            redFound = true;
        }
        if (*iter == Color::Blue) {
            ASSERT(!blueFound);
            blueFound = true;
        }
        if (*iter == Color::Green) {
            ASSERT(!greenFound);
            greenFound = true;
        }
    }
    ASSERT(redFound == hasRed);
    ASSERT(blueFound == hasBlue);
    ASSERT(greenFound == hasGreen);
}

int main() {
    try {
        // Test the getAllValues() method.
        
        ASSERT(Color::getAllValues().size() == 3);
        ASSERT(Color::getAllValues()[Color::Red.getIndex()] == Color::Red);
        ASSERT(Color::getAllValues()[Color::Green.getIndex()] == Color::Green);
        ASSERT(Color::getAllValues()[Color::Blue.getIndex()] == Color::Blue);

        // Verify that acolor got assigned correctly.
        
        ASSERT(acolor == Color::Green);
        ASSERT(acolor != Color::Red);

        // Test creating EnumSets.
        
        EnumSet<Color> set1;
        verifyContents(set1, false, false, false);
        set1 += Color::Red;
        set1 += Color::Blue;
        verifyContents(set1, true, true, false);
        EnumSet<Color> set2;
        verifyContents(set2, false, false, false);
        set2 += Color::Blue;
        verifyContents(set2, false, true, false);
        
        // Test containsAll() and containsAny().
        
        ASSERT(set1.containsAll(set2));
        ASSERT(!set2.containsAll(set1));
        ASSERT(set1.containsAny(set2));
        ASSERT(set2.containsAny(set1));
        
        // Test comparisons of set equality.
        
        ASSERT(set1 != set2);
        set1 -= Color::Red;
        ASSERT(set1.size() == 1);
        ASSERT(set1 == set2);
        
        // Now build two sets we can use for testing lots of operators.
        
        set1.clear();
        set2.clear();
        verifyContents(set1, false, false, false);
        verifyContents(set2, false, false, false);
        set1 += Color::Red;
        set1 += Color::Blue;
        set2 += Color::Blue;
        set2 += Color::Green;
        
        // Test operators.
        
        verifyContents(set1+set2, true, true, true);
        verifyContents(set1|set2, true, true, true);
        verifyContents(set1-set2, true, false, false);
        verifyContents(set2-set1, false, false, true);
        verifyContents(set1&set2, false, true, false);
        verifyContents(set1^set2, true, false, true);
        verifyContents(~set1, false, false, true);
        verifyContents(~set2, true, false, false);
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
