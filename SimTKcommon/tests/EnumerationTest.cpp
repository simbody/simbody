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

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using std::sqrt;
using namespace SimTK;

class Color : public Enumeration<Color> {
public:
    enum Index {RedIndex = 0, GreenIndex = 1, BlueIndex = 2};
    static const Color Red;
    static const Color Green;
    static const Color Blue;
private:
    Color();
    Color(const Color& thisElement, int index, const char* name);
    static void initValues();
    friend class Enumeration<Color>;
};

template class Enumeration<Color>;    // explicit instantiations to make sure we see all the errors
template class EnumerationSet<Color>;

Color acolor = Color::Green; // Verify that accessing a constant that has not yet been initialized works correctly.

const Color Color::Red;
const Color Color::Green;
const Color Color::Blue;

Color::Color() : Enumeration<Color>() {
}

Color::Color(const Color& thisElement, int index, const char* name) : Enumeration<Color>(thisElement, index, name) {
}

void Color::initValues() {
    new(&const_cast<Color&>(Red)) Color(Red, RedIndex, "Red");
    new(&const_cast<Color&>(Green)) Color(Green, GreenIndex, "Green");
    new(&const_cast<Color&>(Blue)) Color(Blue, BlueIndex, "Blue");
}

void verifyContents(EnumerationSet<Color> set, bool hasRed, bool hasBlue, bool hasGreen) {
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
    ASSERT(set.empty() == (size == 0));
    bool redFound = false, blueFound = false, greenFound = false;
    for (EnumerationSet<Color>::iterator iter = set.begin(); iter != set.end(); ++iter) {
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

void verifyName(Enumeration<Color> value) {
    ASSERT(value.getName() == Enumeration<Color>::getValue(value.getIndex()).getName());
}

int main() {
    try {
        // Test the size() and getValue() methods.
        
        ASSERT(Color::size() == 3);
        ASSERT(Color::getValue(Color::Red.getIndex()) == Color::Red);
        ASSERT(Color::getValue(Color::Green.getIndex()) == Color::Green);
        ASSERT(Color::getValue(Color::Blue.getIndex()) == Color::Blue);
        
        // Test iterating over the set of all possible values.
        
        Color color = Color::Red;
        for (Color::iterator iter = Color::begin(); iter != Color::end(); ++iter) {
            ASSERT(*iter == color);
            if (color.getIndex() < Color::size()-1)
                color++;
        }

        // Verify that acolor got assigned correctly.
        
        ASSERT(acolor == Color::Green);
        ASSERT(acolor != Color::Red);
        
        // Test the increment and decrement operators.
        
        color = Color::Red;
        verifyName(++color);
        verifyName(--color);
        verifyName(color++);
        verifyName(color--);

        // Test creating EnumerationSets.
        
        EnumerationSet<Color> set1;
        verifyContents(set1, false, false, false);
        set1 |= Color::Red;
        set1 |= Color::Blue;
        verifyContents(set1, true, true, false);
        EnumerationSet<Color> set2;
        verifyContents(set2, false, false, false);
        set2 |= Color::Blue;
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
        set1 |= Color::Red;
        set1 |= Color::Blue;
        set2 |= Color::Blue;
        set2 |= Color::Green;
        
        // Test operators.
        
        verifyContents(set1|set2, true, true, true);
        verifyContents(set1-set2, true, false, false);
        verifyContents(set2-set1, false, false, true);
        verifyContents(set1&set2, false, true, false);
        verifyContents(set1^set2, true, false, true);
        verifyContents(~set1, false, false, true);
        verifyContents(~set2, true, false, false);
        Color c = Color::Red;
        verifyContents(c++, true, false, false);
        verifyContents(c, false, false, true);
        verifyContents(++c, false, true, false);
        verifyContents(c--, false, true, false);
        verifyContents(c, false, false, true);
        verifyContents(--c, true, false, false);
        
        // Try binary operators involving individual elements.
        
        verifyContents(Color::Blue|Color::Green, false, true, true);
        verifyContents(Color::Green|set1, true, true, true);
        verifyContents(set1|Color::Green, true, true, true);
        verifyContents(Color::Blue&Color::Green, false, false, false);
        verifyContents(Color::Blue&Color::Blue, false, true, false);
        verifyContents(Color::Blue&set1, false, true, false);
        verifyContents(set1&Color::Blue, false, true, false);
        verifyContents(Color::Blue^Color::Green, false, true, true);
        verifyContents(Color::Blue^Color::Blue, false, false, false);
        verifyContents(Color::Blue^set1, true, false, false);
        verifyContents(set1^Color::Blue, true, false, false);
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
