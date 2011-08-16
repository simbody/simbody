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

/**
 * Given the number of values expected and found in a set of bins, verify that the distribution is correct.
 */

void verifyDistribution(int expected[], int found[], int bins) {
    for (int i = 0; i < bins; ++i) {
        Real dev = sqrt((Real) expected[i]);
        ASSERT(found[i] >= expected[i]-4*dev && found[i] <= expected[i]+4*dev)
    }
}

/**
 * Given a set of Reals, verify that they satisfy a uniform distribution between 0 and 1.
 */

void verifyUniformDistribution(Real min, Real max, Real value[], int length) {
    int expected[10], found[10];
    for (int i = 0; i < 10; ++i) {
        expected[i] = length/10;
        found[i] = 0;
    }
    for (int i = 0; i < length; ++i) {
        ASSERT(value[i] >= min)
        ASSERT(value[i] < max)
        int index = (int) ((value[i]-min)*10/(max-min));
        found[index]++;
    }
    verifyDistribution(expected, found, 10);
}

/**
 * Given a set of ints, verify that they satisfy a uniform distribution between 0 and max.
 */

void verifyUniformDistribution(int min, int max, int value[], int length) {
    int range = max-min;
    int* expected = new int[range];
    int* found = new int[range];
    for (int i = 0; i < range; ++i) {
        expected[i] = length/range;
        found[i] = 0;
    }
    for (int i = 0; i < length; ++i) {
        ASSERT(value[i] >= min)
        ASSERT(value[i] < max)
        found[value[i]-min]++;
    }
    verifyDistribution(expected, found, range);
    delete[] expected;
    delete[] found;
}

/**
 * Given a set of values, verify that they satisfy a Gaussian distribution.
 */

void verifyGaussianDistribution(Real mean, Real stddev, Real value[], int length) {
    int expected[6], found[6];
    expected[0] = expected[5] = (int) (0.0228*length);
    expected[1] = expected[4] = (int) (0.1587*length-expected[0]);
    expected[2] = expected[3] = (int) (0.5*length-expected[1]);
    for (int i = 0; i < 6; ++i)
        found[i] = 0;
    for (int i = 0; i < length; ++i) {
        Real val = (value[i]-mean)/stddev;
        if (val < -2)
            found[0]++;
        else if (val < -1)
            found[1]++;
        else if (val < 0)
            found[2]++;
        else if (val < 1)
            found[3]++;
        else if (val < 2)
            found[4]++;
        else
            found[5]++;
    }
    verifyDistribution(expected, found, 6);
}

void testUniform() {
    Random::Uniform rand;
    ASSERT(rand.getMin() == 0.0)
    ASSERT(rand.getMax() == 1.0)

    // Try generating a bunch of random numbers, and make sure they are distributed uniformly between 0 and 1.
    
    Real value[2001];
    value[2000] = 123.4;
    rand.setSeed(1);
    for (int i = 0; i < 2000; ++i)
        value[i] = rand.getValue();
    verifyUniformDistribution(0.0, 1.0, value, 2000);
    
    // Reset the random number generator, and make sure it produces the same values again.
    
    rand.setSeed(1);
    for (int i = 0; i < 2000; ++i)
        ASSERT(value[i] == rand.getValue())
    
    // Now try asking for a whole array at a time, and verify that it still gives the same results.
    
    Real value2[2001];
    value2[2000] = 567.8;
    rand.setSeed(1);
    rand.fillArray(value2, 2000);
    for (int i = 0; i < 2000; ++i)
        ASSERT(value[i] == value2[i])
    
    // Set the seed to a different value, and verify that the results are different.
    
    rand.setSeed(2);
    rand.fillArray(value2, 2000);
    for (int i = 0; i < 2000; ++i)
        ASSERT(value[i] != value2[i])
    
    // Change the range and test the distribution.
    
    rand.setMin(5.0);
    rand.setMax(20.0);
    ASSERT(rand.getMin() == 5.0)
    ASSERT(rand.getMax() == 20.0)
    rand.fillArray(value2, 2000);
    verifyUniformDistribution(5.0, 20.0, value2, 2000);
    
    // Try generating uniform integers.
    
    int value3[2001];
    value3[2000] = -99;
    rand.setSeed(3);
    for (int i = 0; i < 2000; ++i)
        value3[i] = rand.getIntValue();
    verifyUniformDistribution(5, 20, value3, 2000);

    // Verify that if we do not explicitly set the seed, every Random object is initialized with a different seed.
    
    Random::Uniform rand1, rand2;
    rand1.fillArray(value, 2000);
    rand2.fillArray(value2, 2000);
    for (int i = 0; i < 2000; ++i)
        ASSERT(value[i] != value2[i])

    // Make sure none of the above operations has overwritten the final array element.
    
    ASSERT(value[2000] == 123.4)
    ASSERT(value2[2000] = 567.8)
    ASSERT(value3[2000] = -99)
}

void testGaussian() {
    Random::Gaussian rand;
    ASSERT(rand.getMean() == 0.0)
    ASSERT(rand.getStdDev() == 1.0)
    
    // Try generating a bunch of Gaussian random numbers, and check the distribution.
    
    Real value[2001];
    value[2000] = 123.4;
    rand.setSeed(1);
    for (int i = 0; i < 2000; ++i)
        value[i] = rand.getValue();
    verifyGaussianDistribution(0.0, 1.0, value, 2000);
    
    // Try getting a whole array at a time.
    
    Real value2[2001];
    value2[2000] = 567.8;
    rand.setSeed(1);
    rand.fillArray(value2, 2000);
    for (int i = 0; i < 2000; ++i)
        ASSERT(value[i] == value2[i])
    
    // Change the parameters and test the distribution.
    
    rand.setMean(10.0);
    rand.setStdDev(7.0);
    ASSERT(rand.getMean() == 10.0)
    ASSERT(rand.getStdDev() == 7.0)
    rand.fillArray(value2, 2000);
    verifyGaussianDistribution(10.0, 7.0, value2, 2000);
    
    // Make sure none of the above operations has overwritten the final array element.
    
    ASSERT(value[2000] == 123.4)
    ASSERT(value2[2000] = 567.8)
}

int main() {
    try {
        testUniform();
        testGaussian();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


